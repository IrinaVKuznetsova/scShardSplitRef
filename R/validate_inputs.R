#' Coerce a coordinate column to integer, safely
#'
#' Genomic coordinates must always be whole-number `integer`, never
#' `double`. Plain `as.integer()` is not safe to use directly on trusted
#' input because it silently returns `NA` (with only a warning) for values
#' outside R's 32-bit integer range, and silently truncates fractional
#' values. This helper fails loudly instead in both cases so a type problem
#' is caught immediately at the point of reading, rather than surfacing much
#' later as a scientific-notation string or a mysteriously missing row.
#'
#' @param x Numeric vector (integer, double, or character coercible to
#'   numeric) of coordinate values.
#' @param label Character scalar used in error messages to identify which
#'   column failed (e.g. `"RegionStart"`).
#'
#' @returns Integer vector, same length as `x`.
#'
#' @examples
#' \dontrun{
#'   .coerce_integer_coord(c(0, 100000000), "RegionStart")
#'   .coerce_integer_coord(3.5, "RegionStart") # errors: not whole numbers
#' }
#'
#' @dev
.coerce_integer_coord <- function(x, label = "coordinate") {
  if (is.integer(x)) {
    return(x)
  }

  if (!is.numeric(x)) {
    x <- suppressWarnings(as.numeric(x))
  }

  non_integer <- x != trunc(x)
  if (any(non_integer, na.rm = TRUE)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "{.arg {label}} must contain whole-number coordinates.",
        x = "Found {sum(non_integer, na.rm = TRUE)} non-integer value{?s}."
      )
    )
  }

  out <- suppressWarnings(as.integer(x))
  overflowed <- is.na(out) & !is.na(x)

  if (any(overflowed)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "{.arg {label}} contains {sum(overflowed)} value{?s} outside R's \\
         32-bit integer range (+/-{.val {.Machine$integer.max}}).",
        i = "This usually means a coordinate is larger than a single \\
             chromosome/scaffold should be \u2014 check the input file."
      )
    )
  }

  out
}

#' Validate BED and GTF files
#'
#' Reads and validates BED and GTF files, checking for:
#' - File existence
#' - TAB-delimited format with correct number of columns
#' - Valid BED coordinate ranges (start < end)
#'
#' @inheritParams determine_split_regions
#'
#' @returns List with named elements:
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
#'   regions <- system.file("extdata/AlgorithmToy.bed",
#'                    package = "scShardSplitRef",
#'                    mustWork = TRUE)
#'
#'   genes <- system.file("extdata/AlgorithmToy.gtf",
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
#' @dev

.validate_inputs <- function(bed, gtf) {
  .check_files_exist(bed, gtf)

  # .read_and_validate_bed()/.read_and_validate_gtf() validate and parse
  # each file from a single readLines() pass (see their docs below) --
  # this replaces a separate .check_tab_file() structural-validation pass
  # followed by a second, independent read.table(file = path, ...) call,
  # which read every file from disk twice.
  bed <- .read_and_validate_bed(bed)
  gtf <- .read_and_validate_gtf(gtf)

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
#' @returns Invisible NULL. Raises an error if validation fails.
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
  .validate_tab_file_lines(path, min_fields, label, max_fields)
  invisible(TRUE)
}

#' Validate a TAB-delimited file and return its data lines
#'
#' Does the same structural validation as [.check_tab_file()] (extension
#' check, comment/blank-line filtering, field-count range check), but
#' returns the validated, filtered lines instead of `invisible(TRUE)`.
#'
#' This exists so the combined read-and-validate functions
#' ([.read_and_validate_bed()], [.read_and_validate_gtf()]) can parse the
#' file from these already-read, already-filtered lines via
#' `read.table(text = ...)` instead of reading the file from disk a second
#' time. [.check_tab_file()] itself is kept as a separate, unchanged
#' function (rather than refactored to call this and expose its result)
#' so any existing direct callers/tests relying on its exact signature and
#' `invisible(TRUE)` return are unaffected.
#'
#' @inheritParams .check_tab_file
#'
#' @returns Character vector of validated data lines (comments and blank
#'   lines removed). Empty character vector if the file has no data lines.
#'
#' @dev
.validate_tab_file_lines <- function(
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
    return(lines)
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
        x = "Problematic lines: {invalid_idx}.",
        i = "Example: {example_bad}",
        i = "Expected: {expected}"
      ),
      class = "invalid_tab_file_error"
    )
  }

  lines
}

#' Find lines with invalid field counts
#'
#' Identifies line indices where the field count violates min/max constraints.
#'
#' @param field_counts Integer vector: number of fields per line.
#' @param min_fields Integer: minimum required fields.
#' @param max_fields Integer or NA: maximum allowed fields.
#'
#' @returns Integer vector: indices of lines with invalid field counts,
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
#' @returns Character: description such as "3-9 TAB-separated fields" or
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
#' @returns Invisible TRUE if both files exist. Raises error otherwise.
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
#' Kept unchanged, with its original path-based signature, for any direct
#' callers/tests. [.validate_inputs()] no longer calls this directly --
#' it uses [.read_and_validate_gtf()] instead, which validates and parses
#' the file from a single `readLines()` pass rather than reading the file
#' from disk once here and again in [.check_tab_file()].
#'
#' @param path Character: path to GTF file.
#'
#' @returns Data frame with 9 columns named: Chr, source, feature, start,
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

  gtf$start <- .coerce_integer_coord(gtf$start, "start")
  gtf$end <- .coerce_integer_coord(gtf$end, "end")

  return(gtf)
}

#' Read and format BED file
#'
#' Reads a BED file and assigns standard column names to the first three
#' columns. Assumes TAB-delimited format with no header.
#'
#' @param path Character: path to BED file.
#'
#' @returns Data frame with first 3 columns named Chr, RegionStart, RegionEnd,
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
  bed$RegionStart <- .coerce_integer_coord(bed$RegionStart, "RegionStart")
  bed$RegionEnd <- .coerce_integer_coord(bed$RegionEnd, "RegionEnd")
  return(bed)
}

#' Read, validate, and parse a BED file in a single pass
#'
#' Combines the structural validation performed by [.check_tab_file()] with
#' the parsing performed by [.read_and_format_bed()], reading the file from
#' disk exactly once (via [.validate_tab_file_lines()]) and parsing the
#' already-read, already-filtered lines with `read.table(text = ...)`,
#' instead of reading the whole file from disk a second time.
#'
#' @param path Character: path to BED file.
#'
#' @returns Data frame with first 3 columns named Chr, RegionStart, RegionEnd
#'   (integer), plus any additional columns from the input file.
#'
#' @examples
#' \dontrun{
#'   bed_df <- .read_and_validate_bed("regions.bed")
#' }
#'
#' @dev
.read_and_validate_bed <- function(path) {
  lines <- .validate_tab_file_lines(path, min_fields = 3L, label = "BED")

  bed <- utils::read.table(
    text = lines,
    sep = "\t",
    header = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  colnames(bed)[seq_len(3L)] <- c("Chr", "RegionStart", "RegionEnd")
  bed$RegionStart <- .coerce_integer_coord(bed$RegionStart, "RegionStart")
  bed$RegionEnd <- .coerce_integer_coord(bed$RegionEnd, "RegionEnd")
  bed
}

#' Read, validate, and parse a GTF file in a single pass
#'
#' Combines the structural validation performed by [.check_tab_file()] with
#' the parsing performed by [.read_and_format_gtf()], reading the file from
#' disk exactly once (via [.validate_tab_file_lines()]) and parsing the
#' already-read, already-filtered lines with `read.table(text = ...)`,
#' instead of reading the whole file from disk a second time.
#'
#' @param path Character: path to GTF file.
#'
#' @returns Data frame with 9 columns named: Chr, source, feature, start,
#'   end, score, strand, frame, attribute (start/end coerced to integer).
#'
#' @examples
#' \dontrun{
#'   gtf_df <- .read_and_validate_gtf("annotations.gtf")
#' }
#'
#' @dev
.read_and_validate_gtf <- function(path) {
  lines <- .validate_tab_file_lines(path, min_fields = 9L, label = "GTF")

  gtf <- utils::read.table(
    text = lines,
    sep = "\t",
    header = FALSE,
    quote = "",
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
  gtf$start <- .coerce_integer_coord(gtf$start, "start")
  gtf$end <- .coerce_integer_coord(gtf$end, "end")
  gtf
}

#' Validate BED coordinate ranges
#'
#' Checks that all BED regions satisfy the standard requirement that
#' RegionStart < RegionEnd (0-based half-open intervals).
#'
#' @param bed Data frame with RegionStart and RegionEnd columns.
#'
#' @returns Invisible `TRUE` if all ranges are valid. Raises error with
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
        i = "Problematic rows: {invalid_ranges}"
      )
    )
  }

  return(invisible(TRUE))
}
