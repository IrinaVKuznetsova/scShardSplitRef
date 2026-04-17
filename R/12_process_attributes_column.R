#' Filter and write the GTF attribute column
#'
#' Keeps only the requested keys from the 9th GTF column (`attribute`) and
#' writes a final tab-delimited GTF file.
#'
#' @param result A data frame from `process_gtf()` containing standard GTF
#'   columns including `attribute`.
#' @param keep_attributes Character vector of attribute keys to keep, *e.g.*,
#'   `c("gene_id", "transcript_id", "gene_name")`.
#' @param genome_name Genome identifier used in the output filename.
#' @param genome_version Genome version string used in the output filename.
#' @param out_path Directory where the output GTF will be written. Defaults to
#'   the current working directory (`"."`).
#'
#' @returns The function writes the filtered GTF file to disk
#'   at `out_path`.
#'
#' @examples
#' # keep_attributes <- c("gene_id", "transcript_id", "gene_name")
#' # final_gtf_df <- process_attributes_column(
#' #   result = processed_gtf,
#' #   keep_attributes = keep_attributes,
#' #   genome_name = "Rye_Lo7",
#' #   genome_version = "v1p1p1"
#' # )
#' @autoglobal
#' @export

process_attributes_column <- function(
  result,
  keep_attributes,
  genome_name,
  genome_version,
  out_path = "."
) {
  if (!dir.exists(out_path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "Output directory does not exist: {.file {out_path}}"
    )
  }

  filtered <- .filter_attributes(result$attribute, keep_attributes)
  .check_filtered_completeness(filtered, result$attribute)

  result$attribute <- filtered
  output_file <- .build_output_filename(out_path, genome_name, genome_version)

  .write_filtered_gtf(result, output_file)

  cli::cli_inform("{.bold Finished writing processed GTF: {output_file}}")

  invisible(NULL)
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
  filtered <- character(length(attributes))

  for (i in seq_along(attributes)) {
    attr_str <- attributes[i]
    pairs <- strsplit(attr_str, ";\\s*")[[1L]]
    pairs <- pairs[nzchar(pairs)]

    kept_pairs <- character(0)

    for (keep_attr in keep_attributes) {
      pattern <- paste0("^", keep_attr, "\\s+")
      matching <- pairs[grepl(pattern, pairs)]

      if (length(matching) > 0L) {
        kept_pairs <- c(kept_pairs, matching)
      }
    }

    filtered[i] <- paste(kept_pairs, collapse = "; ")
  }

  filtered
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
  empty_idx <- which(!nzchar(filtered) & nzchar(original) > 0L)

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
