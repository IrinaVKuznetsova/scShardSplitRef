#' Process a GTF file by assigning split-region sequence names
#'
#' Reads a GTF and a BED-like "split regions" file, assigns each GTF feature to
#' exactly one split-region interval, and replaces the GTF `Chr` field with a
#' split-region name of the form `Chr:RegionStart-RegionEnd`.
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
#' @param gtf Character: File path to the GTF file.
#' @param split_regions_bed Character: File path to a BED-like file
#'   containing split intervals as produced by [determine_split_regions()]. 
#'   The first three columns must be: `Chr`, `RegionStart`, `RegionEnd`.
#'
#' @return A data frame containing the updated GTF rows, with columns:
#' `Chr`, `source`, `feature`, `start`, `end`, `score`, `strand`, `frame`,
#' `attribute`. The output row order matches the input GTF order.
#'
#' @examples
#' # Minimal working example using temporary files
#'
#'   gtf_file <- system.file("extdata",
#'              "A3_toy_all_scenarios_2chr.gtf",
#'              package = "scShardSplitRef",
#'              mustWork = TRUE)
#'   reg_file <- system.file("extdata",
#'                    "IN0_toy_centromeres_for_gtf.bed",
#'                    package = "scShardSplitRef",
#'                    mustWork = TRUE)
#'
#' out <- process_gtf(gtf, reg_file)
#' out
#'
#' @autoglobal
#' @export
process_gtf <- function(gtf, split_regions_bed) {
  cli::cli_alert_info("Validating inputs and reading files...")
  validated <- validate_inputs(
    gtf_path = gtf,
    bed_path = split_regions_bed
  )
  gtf <- validated$gtf
  regions <- validated$bed

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

  cli::cli_alert_success("Processing complete. Returning formatted GTF.")
  out
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
    "%s:%s-%s",
    regions$Chr,
    regions$RegionStart,
    regions$RegionEnd
  )
  regions
}

#' Identify rows where a feature lies fully within a split region
#'
#' Applies the matching rule used by [process_gtf()].
#'
#' @param df Data frame containing `start`, `end`, `RegionStart`, `RegionEnd`.
#'
#' @return Logical vector of length `nrow(df)`.
#' @dev
.rows_within_region <- function(df) {
  df$start >= df$RegionStart & df$start < df$RegionEnd & df$end <= df$RegionEnd
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
  # no matches is okay here; missing-features check handles that
  if (nrow(matched) == 0L) {
    return(invisible(TRUE))
  }

  dup_ids <- unique(matched$unique_ID[duplicated(matched$unique_ID)])
  n_dup <- length(dup_ids)

  if (n_dup > 0L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "Ambiguous split-region matches detected for GTF feature{?s}.",
        x = "{n_dup} feature ID{?s} matched multiple intervals.",
        i = "Feature ID{?s}: {dup_ids}",
        i = "Ensure the split-regions BED has no overlapping intervals per
             chromosome."
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
  preview <- paste0(
    missing_rows$Chr[seq_len(preview_n)],
    ":",
    missing_rows$feature[seq_len(preview_n)],
    ":",
    missing_rows$start[seq_len(preview_n)],
    ":",
    missing_rows$end[seq_len(preview_n)]
  )

  cli::cli_abort(
    call = rlang::caller_env(),
    c(
      "Some GTF feature{?s} do not fit within any split-region interval.",
      x = "{n_missing} unmatched feature{?s} found.",
      i = "Example{?s}: {preview}",
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
