#' @autoglobal
#' @export

process_attributes_column <- function(
  result,
  keep_attributes,
  out_path,
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

  result
}
