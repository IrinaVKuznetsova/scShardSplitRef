#' Validate GTF and BED files
#'
#' Reads and validates GTF and BED files, checking for:
#' - File existence
#' - TAB-delimited format with correct number of columns
#' - Valid BED coordinate ranges (start < end)
#'
#' @inheritParams determine_split_regions
#'
#' @return List with named elements:
#'   - `bed`: Data frame with at least 3 BED columns plus any
#'      additional columns from the input file.
#'   - `gtf`: Data frame with 9 GTF columns
#'
#' @details
#' Performs four checks in sequence:
#' 1. Both files exist
#' 2. GTF has >= 9 columns
#' 3. BED has >= 3 columns
#' 4. All BED ranges satisfy RegionStart < RegionEnd
#'
#' @examples
#'   regions <- system.file("extdata",
#'                    "IN0_toy_centromeres_for_gtf.bed",
#'                    package = "scShardSplitRef",
#'                    mustWork = TRUE)
#'   genes <- system.file("extdata",
#'              "A3_toy_all_scenarios_2chr.gtf",
#'              package = "scShardSplitRef",
#'              mustWork = TRUE)
#'
#'   result <- validate_inputs(regions, genes)
#'
#'   bed <- result$bed
#'   gtf <- result$gtf
#'
#'   gtf
#'   bed
#'
#' @export
validate_inputs <- function(bed, gtf) {
  .check_files_exist(bed, gtf)

  # check ends with to use .bed or .gtf
  #
  # check end of line encoding
  # check that files are tab delimited

  # mapply() is used here rather than Map since we just need a check, no return list
  mapply(
    .check_tab_file,
    file = list(bed, gtf),
    min_fields = c(3L, 9L),
    label = c("BED", "GTF"),
    SIMPLIFY = FALSE
  )

  bed <- .read_and_format_bed(bed)
  gtf <- .read_and_format_gtf(gtf)

  .validate_bed_ranges(bed)

  cli::cli_inform("{.bold All files loaded and validated successfully.}")

  return(list(bed = bed, gtf = gtf))
}

#' Check TAB-delimited file structure
#'
#' Validates that a file has the correct number of TAB-separated fields.
#' Comment lines starting with '#' are skipped.
#'
#' @param path Character: file path to validate.
#' @param min_fields Integer: minimum number of required fields per line.
#' @param label Character: human-readable file type label (*e.g.*, "GTF", "BED")
#'   for error messages.
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
#'     label = "GTF"
#'   )
#' }
#' @returns Invisible TRUE if valid, else errors if invalid.
#' @dev
.check_tab_file <- function(
  path,
  min_fields,
  label,
  max_fields = NA_integer_
) {
  if (!endsWith(tolower(path), sprintf(".%s", tolower(label)))) {
    cli::cli_abort(c(
      "File extension does not match expected format for {label}.",
      "i" = "Expected extension: {sprintf('.%s', tolower(label))}",
      "x" = "Got: {path}"
    ))
  }

  lines <- readLines(path)
  lines <- lines[!startsWith(lines, "#") & nzchar(trimws(lines))]
  if (length(lines) == 0L) {
    return(invisible(TRUE))
  }

  field_counts <- lengths(stringi::stri_split_fixed(lines, "\t"))

  if (!all(field_counts >= 1L)) {
    cli::cli_abort("File does not appear to be tab-separated.")
  }

  invalid_idx <- .find_invalid_field_counts(
    field_counts,
    min_fields,
    max_fields
  )

  if (length(invalid_idx) > 0L) {
    example_bad <- lines[invalid_idx[1L]]
    expected <- .format_field_expectation(min_fields, max_fields)
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "Invalid {label} file format.",
        x = "Problematic line{?s}: {invalid_idx}.",
        i = "Example: {example_bad}",
        i = "Expected: {expected}"
      )
    )
  }

  return(invisible(TRUE))
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
#'   .format_field_expectation(9L, NA_integer_)
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

#' Check that files exist
#'
#' Validates that both BED and GTF file paths point to existing files.
#'
#' @param bed_path Character: BED file path.
#' @param gtf_path Character: GTF file path.
#'
#' @return Invisible TRUE if both files exist. Raises error otherwise.
#'
#' @examples
#' \dontrun{
#'   .check_files_exist("regions.bed", "genes.gtf")
#' }
#'
#' @dev
.check_files_exist <- function(bed_path, gtf_path) {
  if (!file.exists(bed_path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "BED file does not exist: {.file {bed_path}}"
    )
  }
  if (!file.exists(gtf_path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "GTF file does not exist: {.file {gtf_path}}"
    )
  }
  return(invisible(TRUE))
}

#' Read and format GTF file
#'
#' Reads a GTF file and assigns standard column names. Comment lines are
#' skipped. Assumes TAB-delimited format with no header.
#'
#' @param path Character: path to GTF file.
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
.read_and_format_gtf <- function(path) {
  gtf <- utils::read.table(
    file = path,
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
#' @param path Character: path to BED file.
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
.read_and_format_bed <- function(path) {
  bed <- utils::read.table(
    file = path,
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
      call = rlang::caller_env(),
      c(
        "Invalid BED coordinates.",
        x = "Each row must satisfy {.field RegionStart < RegionEnd}.",
        i = "Problematic row{?s}: {invalid_ranges}"
      )
    )
  }

  return(invisible(TRUE))
}
