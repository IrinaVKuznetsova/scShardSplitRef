#' Filter and write the GTF attribute column
#'
#' Keeps only the requested keys from the 9th GTF column (`attribute`) and
#' writes a final tab-delimited GTF file.
#'
#' @param result A data frame from `process_gtf()` containing standard GTF
#'   columns including `attribute`.
#' @param keep_attributes Character vector of attribute keys to keep, e.g.
#'   `c("gene_id", "transcript_id", "gene_name")`.
#' @param out_path Directory where the output GTF will be written. Defaults to
#'   the current working directory (`"."`).
#' @param genome_name Genome identifier used in the output filename.
#' @param genome_version Genome version string used in the output filename.
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
  out_path = ".",
  genome_name,
  genome_version
) {
  # Regular expression for matching attribute keys
  pat <- sprintf("^(%s)\\s", paste(keep_attributes, collapse = "|"))

  # Vectorised: for each attribute, split, filter by key, glue with ";"
  filtered <- vapply(
    strsplit(result$attribute, ";\\s*"),
    function(parts) {
      subparts <- grep(pat, parts, value = TRUE)
      if (length(subparts)) {
        sprintf("%s;", paste(subparts, collapse = "; "))
      } else {
        ""
      }
    },
    character(1L)
  )

  if (any(filtered == "")) {
    cli::cli_warn(c(
      "Not all attribute fields were found in some rows.",
      "i" = "{sum(filtered == '')} of {length(filtered)} rows are empty after
      filtering."
    ))
  }

  result$attribute <- filtered

  output_file <- file.path(
    out_path,
    sprintf("B1_FINAL_MODIFIED_GTF_%s_%s.gtf", genome_name, genome_version)
  )

  utils::write.table(
    result,
    file = output_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  cli::cli_alert_success("Finished writing processed GTF: {output_file}")

}
