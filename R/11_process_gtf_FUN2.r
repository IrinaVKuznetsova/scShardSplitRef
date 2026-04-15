#' Process the GTF file
#'
#' @param gtf_path file path to the GTF file.
#' @param split_regions_path file path to a BED-like file containing split
#' intervals with chromosome, start, and end columns.
#'
#' @examples
#' #TODO: Add examples of use for this fn.
#'
#' @returns A data frame containing the updated GTF rows.
#' @autoglobal

process_gtf <- function(
  gtf_path,
  split_regions_path
) {
  cli::cli_alert_info("Validating inputs and reading files...")
  validated <- validate_inputs(
    gtf_path = gtf_path,
    bed_path = split_regions_path
  )
  gtf_file <- validated$gtf
  split_regions <- validated$bed

  gtf_file$unique_ID <- sprintf("ID%i", seq_len(nrow(gtf_file)))
  split_regions$NEW <- sprintf(
    "%s:%s-%s",
    split_regions$Chr,
    split_regions$RegionStart,
    split_regions$RegionEnd
  )

  df_merged_chr <- merge(
    gtf_file,
    split_regions,
    by = "Chr",
    all.x = TRUE,
    sort = FALSE
  )

  matched_rows <- with(
    df_merged_chr,
    start >= RegionStart & start < RegionEnd & end <= RegionEnd
  )
  matched <- df_merged_chr[matched_rows, , drop = FALSE]

  duplicate_matches <- table(matched$unique_ID)
  duplicate_matches <- names(duplicate_matches[duplicate_matches > 1L])
  if (length(duplicate_matches) > 0L) {
    cli::cli_abort(c(
      "Ambiguous split-region matches detected for GTF features.",
      "i" = "Feature IDs with multiple matches: {duplicate_matches}",
      "i" = "Ensure the split-regions BED does not contain overlapping
             intervals for the same chromosome."
    ))
  }

  missing_ids <- setdiff(gtf_file$unique_ID, matched$unique_ID)
  if (length(missing_ids) > 0L) {
    missing_rows <- gtf_file[
      gtf_file$unique_ID %in% missing_ids,
      c("Chr", "feature", "start", "end"),
      drop = FALSE
    ]
    missing_preview <- apply(
      utils::head(missing_rows, 5L),
      1L,
      function(row) paste(row, collapse = ":")
    )

    cli::cli_abort(c(
      "Some GTF features do not fit within any split-region interval.",
      "i" = paste("Unmatched features:", length(missing_ids)),
      "i" = paste("Examples:", paste(missing_preview, collapse = "; ")),
      "i" = paste(
        "Check that the split-regions BED fully covers the annotation",
        "coordinates and that boundary-crossing features were handled",
        "before running process_gtf()."
      )
    ))
  }

  cli::cli_alert_info(
    "Assigning split-region sequence names to GTF features..."
  )
  matched <- matched[
    match(gtf_file$unique_ID, matched$unique_ID),
    ,
    drop = FALSE
  ]

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

  cli::cli_alert_success("Processing complete. Returning formatted GTF.")
  out
}

