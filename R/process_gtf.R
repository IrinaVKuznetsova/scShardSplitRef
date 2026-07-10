#' Process a GTF file by assigning split-region sequence names
#'
#' Reads a BED-like "split regions" file and a GTF, assigns each GTF feature to
#' exactly one split-region interval, and replaces the GTF `Chr` field with a
#' split-region name of the form `Chr-RegionStart-RegionEnd`.
#'
#' A feature is considered inside a region when:
#' - `start >= RegionStart`
#' - `start <  RegionEnd`
#' - `end   <= RegionEnd`
#'
#' The function aborts if:
#' - any GTF feature matches more than one region (overlapping regions), or
#' - any GTF feature matches no region (regions don't fully cover the GTF).
#'
#' @section Output file:
#' The function will automatically name the file based on the input file and
#' prepends the following string onto the filename:
#' "B1_FINAL_MODIFIED_GTF_*genome_name*_*genome_version*.gtf", placing
#' the file in the directory specified in the `out_path`.
#'
#' @param split_regions_bed Character: File path to a BED-like file
#'   containing split intervals as produced by [determine_split_regions()].
#'   The first three columns must be: `Chr`, `RegionStart`, `RegionEnd`.
#' @inheritParams determine_split_regions
#' @param keep_attributes Character vector of attribute keys to keep (e.g.,
#'   `c("gene_id", "transcript_id", "gene_name")`). If provided, the GTF
#'   attribute column will be filtered to keep only these fields. Defaults to
#'   `NULL` (no filtering).
#' @param genome_name Genome identifier used in the output filename.
#' @param genome_version Genome version string used in the output filename.
#' @param out_path Directory where the output GTF will be written. Defaults to
#'   the current working directory (`"."`).
#'
#' @returns The function writes the filtered GTF file to disk.
#'
#' @examples
#' reg_file <- system.file(
#'   "extdata/AlgorithmToy_DetermineSplitReg.bed",
#'   package = "scShardSplitRef",
#'   mustWork = TRUE
#' )
#' gtf_file <- system.file(
#'   "extdata/AlgorithmToy.gtf",
#'   package = "scShardSplitRef",
#'   mustWork = TRUE
#' )
#'
#' process_gtf(
#'   split_regions_bed = reg_file,
#'   gtf = gtf_file,
#'   genome_name = "Example",
#'   genome_version = "v1.0",
#'   out_path = tempdir()
#' )
#'
#' @autoglobal
#' @export
process_gtf <- function(
  split_regions_bed,
  gtf,
  genome_name,
  genome_version = NULL,
  keep_attributes = NULL,
  out_path = "."
) {
  cli::cli_alert_info("Validating inputs and reading files...")

  # Validate that genome_name and genome_version are provided
  if (missing(genome_name) || is.null(genome_name)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(x = "Parameter {.arg genome_name} is required.")
    )
  }

  if (is.null(genome_version)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(x = "Parameter {.arg genome_version} is required.")
    )
  }
  validated <- .validate_inputs(
    bed = split_regions_bed,
    gtf = gtf
  )
  regions <- validated$bed
  gtf <- validated$gtf

  gtf$unique_ID <- sprintf("ID%i", seq_len(nrow(gtf)))
  regions <- .prepare_split_regions(regions)

  merged <- merge(
    gtf,
    regions[, c("Chr", "RegionStart", "RegionEnd", "NEW"), drop = FALSE],
    by = "Chr",
    all.x = TRUE,
    sort = FALSE
  )

  matched <- merged[.rows_within_region(merged), , drop = FALSE]

  .abort_if_ambiguous_matches(matched)
  .abort_if_unmatched_features(matched, gtf)

  matched$start <- matched$start - matched$RegionStart + 1
  matched$end <- matched$end - matched$RegionStart + 1

  matched <- .restore_gtf_order(matched, gtf$unique_ID)

  out <- .build_gtf_output(matched)

  # Handle attribute filtering if parameters are provided
  if (!is.null(keep_attributes)) {
    filtered <- .filter_attributes(out$attribute, keep_attributes)
    .check_filtered_completeness(filtered, out$attribute)
    out$attribute <- filtered
  }

  # Write to file
  if (!dir.exists(out_path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "Output directory does not exist: {.file {out_path}}"
    )
  }
  output_file <- .build_output_filename(out_path, genome_name, genome_version)
  .write_filtered_gtf(out, output_file)
  cli::cli_alert_info(c(
    i = "{.bold Finished writing processed GTF: {output_file}}"
  ))

  invisible(NULL)
}

#' Prepare split regions by adding `NEW` region names
#'
#' Adds a `NEW` column to the split-regions data frame, formatted as
#' `Chr:RegionStart-RegionEnd`.
#'
#' @param regions Data frame with columns `Chr`, `RegionStart`, `RegionEnd`.
#'
#' @return The input data frame with an added `NEW` column.
#' @dev
.prepare_split_regions <- function(regions) {
  regions$NEW <- sprintf(
    "%s-%s-%s",
    regions$Chr,
    regions$RegionStart,
    regions$RegionEnd
  )
  regions
}

#' Identify rows where a feature lies fully within a split region
#'
#' Applies the matching rule used by [process_gtf()]. A feature `start`, `end`
#' should be fully inside a region:
#'  `RegionStart`, `RegionEnd`:
#'  `start >= RegionStart`
#'  `end <= RegionEnd`
#'
#'
#' @param df Data frame containing `start`, `end`, `RegionStart`, `RegionEnd`.
#'
#' @return Logical vector of length `nrow(df)`.
#' @dev
.rows_within_region <- function(df) {
  df$start >= df$RegionStart & df$end <= df$RegionEnd
}

#' Abort if any feature matches multiple split regions
#'
#' Detects duplicated `unique_ID` values in the matched table, which indicates
#' overlapping split regions (or otherwise ambiguous assignments).
#'
#' @param matched Data frame containing a `unique_ID` column.
#'
#' @return Invisibly `TRUE` if ok; otherwise aborts.
#' @dev

.abort_if_ambiguous_matches <- function(matched) {
  if (nrow(matched) == 0L) {
    return(invisible(TRUE))
  }

  dup_ids <- unique(matched$unique_ID[duplicated(matched$unique_ID)])
  n_dup <- length(dup_ids)

  if (n_dup > 0L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "{n_dup} ambiguous split-region matches detected.",
        x = "{n_dup} feature IDs matched multiple intervals.",
        i = "Affected IDs:",
        i = paste0("* {.val ", dup_ids, "}"),
        i = "Ensure the split-regions BED has no overlapping intervals per chromosome."
      ),
      .envir = environment()
    )
  }

  invisible(TRUE)
}

#' Abort if any GTF features match no split region
#'
#' @param matched Data frame containing `unique_ID`.
#' @param gtf Original GTF data frame containing `unique_ID`, `Chr`, `feature`,
#'   `start`, and `end`.
#'
#' @return Invisibly `TRUE` if ok; otherwise aborts.
#' @dev

.abort_if_unmatched_features <- function(matched, gtf) {
  missing_ids <- setdiff(gtf$unique_ID, matched$unique_ID)
  n_missing <- length(missing_ids)

  if (n_missing == 0L) {
    return(invisible(TRUE))
  }

  missing_rows <- gtf[
    gtf$unique_ID %in% missing_ids,
    c("Chr", "feature", "start", "end"),
    drop = FALSE
  ]

  preview_n <- min(5L, nrow(missing_rows))
  preview_items <- sprintf(
    "{.code %s:%s:%s-%s}",
    missing_rows$Chr[seq_len(preview_n)],
    missing_rows$feature[seq_len(preview_n)],
    missing_rows$start[seq_len(preview_n)],
    missing_rows$end[seq_len(preview_n)]
  )

  cli::cli_abort(
    call = rlang::caller_env(),
    c(
      "{n_missing} GTF features do not fall within any split-region interval.",
      x = "{n_missing} unmatched feature{?s} found.",
      i = "Example:",
      i = paste0("* ", preview_items),
      i = "Ensure split regions fully cover the annotation coordinates.",
      i = "Handle boundary-crossing features before calling {.fn process_gtf}."
    ),
    .envir = environment()
  )
}

#' Restore the original GTF row order
#'
#' @param matched Data frame containing `unique_ID`.
#' @param gtf_unique_ids Character vector of `unique_ID` values in the original
#'   GTF order.
#'
#' @return `matched` reordered to match the input GTF order.
#' @dev
.restore_gtf_order <- function(matched, gtf_unique_ids) {
  matched[match(gtf_unique_ids, matched$unique_ID), , drop = FALSE]
}

#' Build the output GTF data frame
#'
#' Selects and renames columns so the output matches GTF layout, replacing
#' `Chr` with the split-region name stored in `NEW`.
#'
#' @param matched Data frame containing `NEW` and standard GTF columns.
#'
#' @return A data frame with GTF columns where the first column is `Chr`.
#' @dev
.build_gtf_output <- function(matched) {
  out <- matched[, c(
    "NEW",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )]

  colnames(out)[1] <- "Chr"
  out
}

#' Filter GTF attributes
#'
#' Extract and keep only specified attributes from GTF attribute column
#'
#' @param attributes Character vector of GTF attribute strings
#' @param keep_attributes Character vector of attribute names to keep
#'
#' @return Character vector of filtered attribute strings
#'
#' @details
#' Extracts key-value pairs from GTF attribute column and keeps only
#' specified attributes in the original order.
#'
#' @dev
.filter_attributes <- function(attributes, keep_attributes) {
  pattern <- paste0("^(?:", paste(keep_attributes, collapse = "|"), ")\\s+")

  str_list <- stringi::stri_split_regex(attributes, ";\\s*")

  vapply(
    str_list,
    FUN.VALUE = character(1L),
    FUN = function(pairs) {
      pairs <- pairs[nzchar(pairs)]
      keep <- stringi::stri_detect_regex(pairs, pattern)
      if (!any(keep)) {
        return("")
      }
      paste(pairs[keep], collapse = "; ")
    }
  )
}

#' Check filtered attributes completeness
#'
#' Verify that filtering didn't result in empty attributes
#'
#' @param filtered Character vector of filtered attributes
#' @param original Character vector of original attributes
#'
#' @return Invisible TRUE if valid
#'
#' @dev
.check_filtered_completeness <- function(filtered, original) {
  is_empty_filtered <- !nzchar(filtered)
  is_nonempty_original <- nzchar(original)

  empty_idx <- which(is_empty_filtered & is_nonempty_original)

  if (length(empty_idx) > 0L) {
    cli::cli_warn(
      c(
        "{length(empty_idx)} lines had no matching attributes.",
        i = "Check your keep_attributes list."
      )
    )
  }

  invisible(TRUE)
}

#' Build output filename for filtered GTF
#'
#' @param out_path Character: output directory
#' @param genome_name Character: genome name
#' @param genome_version Character: genome version
#'
#' @return Character: full file path
#'
#' @dev
.build_output_filename <- function(out_path, genome_name, genome_version) {
  filename <- sprintf(
    "B1_FINAL_MODIFIED_GTF_%s_%s.gtf",
    genome_name,
    genome_version
  )
  file.path(out_path, filename)
}

#' Write filtered GTF file
#'
#' @param result Data frame with filtered GTF data
#' @param output_file Character: output file path
#'
#' @return Called for its side-effects of writing a file, returns an invisible
#'  `TRUE`.
#'
#' @dev
.write_filtered_gtf <- function(result, output_file) {
  invisible(utils::write.table(
    result,
    file = output_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  ))
  invisible(TRUE)
}
