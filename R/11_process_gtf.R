#' Process the GTF file
#'
#' Assigns GTF features to split genomic regions by matching feature
#' coordinates with split interval boundaries. Each feature must fit entirely
#' within exactly one split region interval.
#'
#' @param gtf_path Character: file path to the GTF file.
#' @param split_regions_path Character: file path to BED-like split intervals
#'   with chromosome, start, and end columns.
#'
#' @return Data frame with GTF data updated with split region identifiers
#'   (formatted as "chr:start-end") in place of the original chromosome column.
#'   Columns are: Chr (new region name), source, feature, start, end, score,
#'   strand, frame, attribute.
#'
#' @details
#' The function:
#' 1. Validates input files (via `validate_inputs()`)
#' 2. Assigns unique IDs to GTF features for tracking
#' 3. Formats split region names as "chr:start-"
#' 4. Matches each feature to the region it falls within
#' 5. Validates that matches are unambiguous (no feature matches multiple regions)
#' 6. Validates that all features match exactly one region
#' 7. Returns GTF with updated Chr column containing region names
#'
#' @examples
#' \dontrun{
#'   result <- process_gtf("annotations.gtf", "split_regions.bed")
#'   head(result)
#' }
#'
#' @export
process_gtf <- function(
  gtf_path,
  split_regions_path
) {
  cli::cli_inform("Validating inputs and reading files.../n")

  validated <- validate_inputs(
    gtf_path = gtf_path,
    bed_path = split_regions_path
  )

  gtf_file <- validated$gtf
  split_regions <- validated$bed

  gtf_file$unique_ID <- sprintf("ID%i", seq_len(nrow(gtf_file)))
  split_regions$NEW <- .format_region_name(
    split_regions$Chr,
    split_regions$RegionStart,
    split_regions$RegionEnd
  )

  matched <- .match_features_to_regions(gtf_file, split_regions)
  .check_match_validity(matched, gtf_file)

  cli::cli_inform("Assigning split-region names to GTF features.../n")

  out <- .build_output_gtf(matched, gtf_file)

  cli::cli_inform("{.bold Processing complete. Returning formatted GTF.}")
  out
}

# ============================================================================
# Core Matching Logic
# ============================================================================

#' Format split region name
#'
#' Constructs a standardized split region identifier in the format
#' "chr:start-end".
#'
#' @param chr Character: chromosome name.
#' @param start Integer: region start coordinate.
#' @param end Integer: region end coordinate.
#'
#' @return Character: formatted region name such as "chr1:100000-200000".
#'
#' @examples
#' \dontrun{
#'   name <- .format_region_name("chr1", 100000, 200000)
#'   # Returns: "chr1:100000-200000"
#' }
#'
#' @dev
.format_region_name <- function(chr, start, end) {
  sprintf("%s:%s-%s", chr, start, end)
}

#' Match GTF features to split regions
#'
#' For each GTF feature, finds the split region that contains it. A feature
#' is considered to match a region if its start coordinate >= region start,
#' start < region end, and end <= region end (entirely contained).
#'
#' @param gtf_file Data frame with GTF data including columns chr, start, end.
#' @param split_regions Data frame with split region intervals including
#'   columns Chr, RegionStart, RegionEnd, and NEW (formatted region names).
#'
#' @return Data frame: merged GTF and split region data for all matching
#'   feature-region pairs. Rows are restricted to valid matches where
#'   features fit entirely within regions.
#'
#' @details
#' Performs an inner merge on chromosome, then filters to keep only rows
#' where the GTF feature coordinates satisfy the containment requirement.
#'
#' @examples
#' \dontrun{
#'   matched <- .match_features_to_regions(gtf_file, split_regions)
#' }
#'
#' @dev
.match_features_to_regions <- function(gtf_file, split_regions) {
  # Merge GTF with split regions on chromosome
  df_merged <- merge(
    gtf_file,
    split_regions,
    by = "Chr",
    all.x = TRUE,
    sort = FALSE
  )

  # Find rows where feature coordinates fit within region
  matched_rows <- with(
    df_merged,
    start >= RegionStart & start < RegionEnd & end <= RegionEnd
  )

  df_merged[matched_rows, , drop = FALSE]
}

#' Check match validity and completeness
#'
#' Validates that feature-region matches satisfy two critical conditions:
#' 1. No feature matches multiple regions (ambiguous matches)
#' 2. All GTF features have exactly one match
#'
#' @param matched Data frame with matched feature-region pairs.
#' @param gtf_file Data frame with original GTF data (for reference).
#'
#' @return Invisible TRUE if valid. Raises error otherwise.
#'
#' @examples
#' \dontrun{
#'   .check_match_validity(matched, gtf_file)
#' }
#'
#' @dev
.check_match_validity <- function(matched, gtf_file) {
  .check_no_ambiguous_matches(matched)
  .check_all_features_matched(matched, gtf_file)
  return(invisible(TRUE))
}

#' Check for ambiguous feature-region matches
#'
#' Ensures that each GTF feature matches to at most one split region.
#' Raises an error if any feature has multiple matches.
#'
#' @param matched Data frame with matched pairs, including unique_ID column.
#'
#' @return Invisible TRUE if no ambiguities. Raises error listing duplicate
#'   feature IDs otherwise.
#'
#' @examples
#' \dontrun{
#'   .check_no_ambiguous_matches(matched)
#' }
#'
#' @dev
.check_no_ambiguous_matches <- function(matched) {
  if (nrow(matched) == 0L) {
    return(invisible(TRUE))
  }

  dup_table <- table(matched$unique_ID)
  duplicates <- names(dup_table[dup_table > 1L])

  if (length(duplicates) > 0L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "Ambiguous split-region matches detected.",
        x = "GTF feature IDs with multiple matches: {duplicates}",
        i = "Ensure split-region BED has non-overlapping intervals."
      )
    )
  }

  return(invisible(TRUE))
}

#' Check that all features are matched
#'
#' Validates that every GTF feature has been assigned to exactly one split
#' region. Raises an informative error if any features are unmatched.
#'
#' @param matched Data frame with successfully matched feature-region pairs.
#' @param gtf_file Data frame with original GTF data.
#'
#' @return Invisible TRUE if all features are matched. Raises error with
#'   examples of unmatched features otherwise.
#'
#' @details
#' Identifies features from gtf_file that are missing from matched data,
#' then reports their count and shows up to 5 examples.
#'
#' @examples
#' \dontrun{
#'   .check_all_features_matched(matched, gtf_file)
#' }
#'
#' @dev
.check_all_features_matched <- function(matched, gtf_file) {
  missing_ids <- setdiff(gtf_file$unique_ID, matched$unique_ID)

  if (length(missing_ids) == 0L) {
    return(invisible(TRUE))
  }

  missing_rows <- gtf_file[
    gtf_file$unique_ID %in% missing_ids,
    c("Chr", "feature", "start", "end"),
    drop = FALSE
  ]

  # Use apply() directly - MARGIN=1 for rows
  missing_preview <- apply(
    utils::head(missing_rows, 5L),
    MARGIN = 1L,
    FUN = \(row) paste(row, collapse = ":")
  )
  cli::cli_abort(
    call = rlang::caller_env(),
    c(
      "Some GTF features do not fit within any split-region interval.",
      x = "Unmatched features: {length(missing_ids)}",
      i = "Examples: {toString(missing_preview)}",
      i = "Ensure split-regions BED fully covers annotation coordinates."
    )
  )
  return(invisible(TRUE))
}

#' Build output GTF data frame
#'
#' Constructs the final GTF output by reordering matched data to match the
#' original GTF order and selecting the appropriate columns.
#'
#' @param matched Data frame with matched feature-region pairs.
#' @param gtf_file Data frame with original GTF data (for ordering reference).
#'
#' @return Data frame with 9 GTF columns: Chr (containing new region names),
#'   source, feature, start, end, score, strand, frame, attribute.
#'   Row order matches the original gtf_file.
#'
#' @examples
#' \dontrun{
#'   out <- .build_output_gtf(matched, gtf_file)
#' }
#'
#' @dev
#'

.build_output_gtf <- function(matched, gtf_file) {
  # Ensure NEW column exists
  if (!"NEW" %in% colnames(matched)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "NEW column not found in matched data."
    )
  }

  # Sort by NEW column (region name)
  matched <- matched[order(matched$NEW), ]

  # Create Chr column from NEW
  matched$Chr <- matched$NEW

  # Select and order columns with Chr first
  result <- matched[,
    c(
      "Chr",
      "source",
      "feature",
      "start",
      "end",
      "score",
      "strand",
      "frame",
      "attribute",
      "NEW"
    ),
    drop = FALSE
  ]

  return(result)
}
