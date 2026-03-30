#' Helper: check TAB-delimited structure (outside main function, for reuse)
#'
#' @dev
#' @autoglobal

check_tab_file <- function(
  path,
  min_fields,
  max_fields = NA_integer_,
  label,
  error_num
) {
  lines <- suppressWarnings(readLines(path))
  lines <- lines[!startsWith(lines, "#")]
  if (length(lines) == 0L) {
    return(invisible(NULL))
  }
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  field_counts <- vapply(split_lines, lengths, integer(1))
  bad_min <- which(field_counts < min_fields)
  bad_max <- if (!is.na(max_fields)) {
    which(field_counts > max_fields)
  } else {
    integer(0)
  }
  bad <- sort(unique(c(bad_min, bad_max)))
  if (length(bad) > 0L) {
    expected <- if (!is.na(max_fields)) {
      glue::glue("{min_fields}-{max_fields} TAB-separated fields")
    } else {
      glue::glue("= {min_fields} TAB-separated fields")
    }
    cli::cli_abort(
      c(
        "ERROR {.val {error_num}}: Invalid {label} formatting.",
        "Problematic line numbers: {.val {paste(bad, collapse = ', ')}}",
        "Example problematic line:",
        "  {lines[bad[1]]}",
        "Expected: {expected}",
        "This file must be TAB-delimited."
      ),
      .envir = environment()
    )
  }
  invisible(NULL)
}

#' Helper: Validate GTF and BED files
#'
#' @param gtf_path User provided file path to a GTF file for validation.
#' @param centromere_path User provided file path to a BED file for validation.
#'
#' @returns A list object containing validated GTF and BED files.
#' @dev
#' @autoglobal
validate_inputs <- function(gtf_path, centromere_path) {
  # --- ERROR 0: Check file paths ----------------------------------------------
  if (!file.exists(gtf_path)) {
    cli::cli_abort("ERROR 0: GTF file does not exist. Check: {gtf_path}")
  }
  if (!file.exists(centromere_path)) {
    cli::cli_abort(
      "ERROR 1: Centromere BED file does not exist. Check: {centromere_path}"
    )
  }

  # --- ERROR 2: Validate GTF --------------------------------------------------
  check_tab_file(
    path = gtf_path,
    min_fields = 9,
    max_fields = NA_integer_,
    label = "GTF",
    error_num = 2L
  )
  # --- ERROR 3: Validate BED --------------------------------------------------
  check_tab_file(
    path = centromere_path,
    min_fields = 3L,
    max_fields = 3L,
    label = "BED",
    error_num = 3L
  )

  # --- Read validated files ---------------------------------------------------
  gtf <- utils::read.table(
    gtf_path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE
  )
  colnames(gtf) <- c(
    "Chr",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )

  bed <- utils::read.table(
    centromere_path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(bed) <- c("Chr", "CentrStart", "ChrEND")

  cli::cli_alert_success("All files loaded and validated successfully.")

  return(list(gtf = gtf, bed = bed))
}
