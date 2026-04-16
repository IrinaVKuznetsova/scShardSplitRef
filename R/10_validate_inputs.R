#' Check TAB-delimited file structure
#'
#' Validates that a file has the correct number of TAB-separated fields.
#' Comment lines starting with '#' are skipped.
#'
#' @param path Character: file path to validate.
#' @param min_fields Integer: minimum number of required fields per line.
#' @param label Character: human-readable file type label (*e.g.*, "GTF", "BED")
#'   for error messages.
#' @param error_num Integer: error code number for identification in messages.
#' @param max_fields Integer or `NA`: maximum number of allowed fields.
#'   If `NA` (default), no upper limit is enforced.
#' @return Invisible NULL. Raises an error if validation fails.
#'
#' @details
#' Reads all lines from the file, skips comments, then checks that each
#' line has a field count in the valid range \[min_fields, max_fields\].
#' Reports line numbers and examples of problematic lines.
#'
#' @examples
#' \dontrun{
#'   .check_tab_file(
#'     "annotations.gtf",
#'     min_fields = 9L,
#'     label = "GTF",
#'     error_num = 2L
#'   )
#' }
#'
#' @dev
.check_tab_file <- function(
  path,
  min_fields,
  label,
  error_num,
  max_fields = NA_integer_
) {
  lines <- readLines(path)
  lines <- lines[!startsWith("#", lines)]
  lines <- lines[nzchar(lines)]

  if (length(lines) == 0L) {
    return(invisible(TRUE))
  }

  split_lines <- strsplit(lines, "\t")
  field_counts <- lengths(split_lines)

  invalid_idx <- .find_invalid_field_counts(
    field_counts,
    min_fields,
    max_fields
  )

  if (length(invalid_idx) > 0L) {
    cli::cli_abort(
      "{error_num}: Invalid {label} file format at line{?s} {invalid_idx}"
    )
  }

  invisible(TRUE)
}

#' Find lines with invalid field counts
#'
#' Identifies line indices where the field count violates min/max constraints.
#'
#' @param field_counts Integer vector: number of fields per line.
#' @param min_fields Integer: minimum required fields.
#' @param max_fields Integer or NA: maximum allowed fields.
#'
#' @return Integer vector: indices of lines with invalid field counts,
#'   sorted and deduplicated. Empty vector if all lines are valid.
#'
#' @examples
#' \dontrun{
#'   bad <- .find_invalid_field_counts(c(2, 3, 2), 3L, NA_integer_)
#' }
#'
#' @dev
.find_invalid_field_counts <- function(
  field_counts,
  min_fields,
  max_fields = NA_integer_
) {
  invalid <- which(field_counts < min_fields)

  if (!is.na(max_fields)) {
    invalid <- c(invalid, which(field_counts > max_fields))
  }

  sort(unique(invalid))
}

#' Report field count validation error
#'
#' Generates a detailed error message about field count violations,
#' including examples, expected format, and advice.
#'
#' @param bad_lines Integer vector: line numbers with invalid field counts.
#' @param lines Character vector: all file lines.
#' @param min_fields Integer: minimum required fields.
#' @param max_fields Integer or NA: maximum allowed fields.
#' @param label Character: file type label.
#' @param error_num Integer: error code number.
#'
#' @return Raises an error. Does not return normally.
#'
#' @examples
#' \dontrun{
#'   .report_field_error(c(1, 3), lines, 3L, NA_integer_, "GTF", 2L)
#' }
#'
#' @dev
.report_field_error <- function(
  bad_lines,
  lines,
  min_fields,
  max_fields,
  label,
  error_num
) {
  expected <- .format_field_expectation(min_fields, max_fields)
  example_bad <- lines[bad_lines[1L]]

  cli::cli_abort(
    c(
      "{.val {error_num}}: Invalid {label} formatting.",
      x = "Problematic line numbers: {bad_lines}.",
      i = "Example: {example_bad}",
      i = "Expected: {expected}",
      i = "File must be TAB-delimited."
    ),
    .envir = environment()
  )
}

#' Format field requirement description
#'
#' Creates a human-readable string describing the allowed field count range.
#'
#' @param min_fields Integer: minimum fields.
#' @param max_fields Integer or NA: maximum fields. If NA, no upper limit.
#'
#' @return Character: description such as "3-9 TAB-separated fields" or
#'   "at least 9 TAB-separated fields".
#'
#' @examples
#' \dontrun{
#'   .format_field_expectation(3L, 9L)
#'   # Returns: "3-9 TAB-separated fields"
#'   format_field_expectation(9L, NA_integer_)
#'   # Returns: "at least 9 TAB-separated fields"
#' }
#'
#' @dev
.format_field_expectation <- function(min_fields, max_fields) {
  if (!is.na(max_fields)) {
    sprintf("%d-%d TAB-separated fields", min_fields, max_fields)
  } else {
    sprintf("at least %d TAB-separated fields", min_fields)
  }
}

# ============================================================================

#' Validate GTF and BED files
#'
#' Reads and validates GTF and BED files, checking for:
#' - File existence
#' - TAB-delimited format with correct number of columns
#' - Valid BED coordinate ranges (start < end)
#'
#' @param gtf_path Character: file path to GTF file.
#' @param bed_path Character: file path to BED file.
#'
#' @return List with named elements:
#'   - `gtf`: Data frame with 9 GTF columns (Chr, source, feature, start,
#'     end, score, strand, frame, attribute)
#'   - `bed`: Data frame with at least 3 BED columns (Chr, RegionStart,
#'      RegionEnd) plus any additional columns from the input file.
#'
#' @details
#' Performs four checks in sequence:
#' 1. Both files exist (ERROR 0 for GTF, ERROR 1 for BED)
#' 2. GTF has >= 9 columns (ERROR 2)
#' 3. BED has >= 3 columns (ERROR 3)
#' 4. All BED ranges satisfy RegionStart < RegionEnd (ERROR 4)
#'
#' @examples
#' \dontrun{
#'   result <- .validate_inputs("genes.gtf", "regions.bed")
#'   gtf <- result$gtf
#'   bed <- result$bed
#' }
#'
#' @export
validate_inputs <- function(gtf_path, bed_path) {
  .check_files_exist(gtf_path, bed_path)
  .check_tab_file(gtf_path, min_fields = 9L, label = "GTF", error_num = 2L)
  .check_tab_file(bed_path, min_fields = 3L, label = "BED", error_num = 3L)

  gtf <- .read_and_format_gtf(gtf_path)
  bed <- .read_and_format_bed(bed_path)

  .validate_bed_ranges(bed)

  cli::cli_inform("All files loaded and validated successfully.")

  return(list(gtf = gtf, bed = bed))
}

#' Check that files exist
#'
#' Validates that both GTF and BED file paths point to existing files.
#'
#' @param gtf_path Character: GTF file path.
#' @param bed_path Character: BED file path.
#'
#' @return Invisible TRUE if both files exist. Raises error otherwise.
#'
#' @examples
#' \dontrun{
#'   .check_files_exist("genes.gtf", "regions.bed")
#' }
#'
#' @dev
.check_files_exist <- function(gtf_path, bed_path) {
  if (!file.exists(gtf_path)) {
    cli::cli_abort(
      "GTF file does not exist: {.file {gtf_path}}"
    )
  }
  if (!file.exists(bed_path)) {
    cli::cli_abort(
      "BED file does not exist: {.file {bed_path}}"
    )
  }
  return(invisible(TRUE))
}

#' Read and format GTF file
#'
#' Reads a GTF file and assigns standard column names. Comment lines are
#' skipped. Assumes TAB-delimited format with no header.
#'
#' @param gtf_path Character: path to GTF file.
#'
#' @return Data frame with 9 columns named: Chr, source, feature, start,
#'   end, score, strand, frame, attribute.
#'
#' @examples
#' \dontrun{
#'   gtf_df <- .read_and_format_gtf("annotations.gtf")
#' }
#'
#' @dev
.read_and_format_gtf <- function(gtf_path) {
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

  return(gtf)
}

#' Read and format BED file
#'
#' Reads a BED file and assigns standard column names to the first three
#' columns. Assumes TAB-delimited format with no header.
#'
#' @param bed_path Character: path to BED file.
#'
#' @return Data frame with first 3 columns named Chr, RegionStart, RegionEnd,
#'   plus any additional columns from the input file.
#'
#' @examples
#' \dontrun{
#'   bed_df <- .read_and_format_bed("regions.bed")
#' }
#'
#' @dev
.read_and_format_bed <- function(bed_path) {
  bed <- utils::read.table(
    bed_path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE
  )

  colnames(bed)[seq_len(3L)] <- c("Chr", "RegionStart", "RegionEnd")
  return(bed)
}

#' Validate BED coordinate ranges
#'
#' Checks that all BED regions satisfy the standard requirement that
#' RegionStart < RegionEnd (0-based half-open intervals).
#'
#' @param bed Data frame with RegionStart and RegionEnd columns.
#'
#' @return Invisible TRUE if all ranges are valid. Raises error with
#'   problematic row numbers otherwise.
#'
#' @examples
#' \dontrun{
#'   .validate_bed_ranges(bed_df)
#' }
#'
#' @dev
.validate_bed_ranges <- function(bed) {
  invalid_ranges <- which(bed$RegionStart >= bed$RegionEnd)

  if (length(invalid_ranges) > 0L) {
    cli::cli_abort(
      c(
        "Invalid BED coordinates.",
        x = "Each row must satisfy {.field RegionStart < RegionEnd}.",
        i = "Problematic row{?s}: {invalid_ranges}"
      )
    )
  }

  return(invisible(TRUE))
}
