# Tests for 10_validate_inputs_FUN1.r

# .find_invalid_field_counts tests
test_that(".find_invalid_field_counts identifies violations", {
  # Min violations
  field_counts <- c(2L, 3L, 2L, 4L)
  expect_identical(.find_invalid_field_counts(field_counts, 3L), c(1L, 3L))

  # Max violations
  field_counts <- c(3L, 5L, 3L, 6L)
  expect_identical(.find_invalid_field_counts(field_counts, 2L, 4L), c(2L, 4L))

  # Both
  field_counts <- c(1L, 3L, 5L, 2L)
  expect_identical(.find_invalid_field_counts(field_counts, 2L, 4L), c(1L, 3L))
})

test_that(".find_invalid_field_counts returns empty when valid", {
  field_counts <- c(3L, 4L, 3L, 5L)
  expect_identical(
    .find_invalid_field_counts(field_counts, 2L, 6L),
    integer(0L)
  )
})

test_that(".find_invalid_field_counts handles NA max", {
  field_counts <- c(3L, 100L, 50L)
  expect_identical(
    .find_invalid_field_counts(field_counts, 2L, NA_integer_),
    integer(0L)
  )
})

test_that(".find_invalid_field_counts deduplicates and sorts", {
  field_counts <- c(1L, 1L, 5L)
  expect_identical(
    .find_invalid_field_counts(field_counts, 2L, 4L),
    c(1L, 2L, 3L)
  )
})

# .format_field_expectation tests
test_that(".format_field_expectation formats correctly", {
  expect_identical(
    .format_field_expectation(3L, 9L),
    "3-9 TAB-separated fields"
  )
  expect_identical(
    .format_field_expectation(9L, NA_integer_),
    "at least 9 TAB-separated fields"
  )
})

# .check_files_exist tests
test_that(".check_files_exist validates both files", {
  gtf_tmp <- tempfile(fileext = ".gtf")
  bed_tmp <- tempfile(fileext = ".bed")
  file.create(gtf_tmp)
  file.create(bed_tmp)
  on.exit({
    unlink(gtf_tmp)
    unlink(bed_tmp)
  })

  expect_invisible(.check_files_exist(gtf_tmp, bed_tmp))
})

test_that(".check_files_exist reports missing files", {
  expect_error(
    .check_files_exist("/nonexistent/file.gtf", tempfile()),
    regexp = "ERROR 0"
  )

  gtf_tmp <- tempfile(fileext = ".gtf")
  file.create(gtf_tmp)
  on.exit(unlink(gtf_tmp))

  expect_error(
    .check_files_exist(gtf_tmp, "/nonexistent/file.bed"),
    regexp = "ERROR 1"
  )
})

# .read_and_format_gtf tests
test_that(".read_and_format_gtf reads and formats", {
  tmpfile <- tempfile(fileext = ".gtf")
  writeLines(
    'chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1";',
    tmpfile
  )
  on.exit(unlink(tmpfile))

  result <- .read_and_format_gtf(tmpfile)
  expect_identical(
    colnames(result),
    c(
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
  )
})

test_that(".read_and_format_gtf skips comments", {
  tmpfile <- tempfile(fileext = ".gtf")
  writeLines(
    c("# Comment", 'chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1";'),
    tmpfile
  )
  on.exit(unlink(tmpfile))

  result <- .read_and_format_gtf(tmpfile)
  expect_identical(nrow(result), 1L)
})

# .read_and_format_bed tests
test_that(".read_and_format_bed reads and formats", {
  tmpfile <- tempfile(fileext = ".bed")
  writeLines("chr1\t0\t100\textra", tmpfile)
  on.exit(unlink(tmpfile))

  result <- .read_and_format_bed(tmpfile)
  expect_identical(colnames(result)[1:3], c("Chr", "RegionStart", "RegionEnd"))
  expect_identical(colnames(result)[4], "V4")
})

test_that(".read_and_format_bed preserves extra columns", {
  tmpfile <- tempfile(fileext = ".bed")
  writeLines(
    c("chr1\t0\t100\tname\tscore", "chr2\t100\t200\tname2\tscore2"),
    tmpfile
  )
  on.exit(unlink(tmpfile))

  result <- .read_and_format_bed(tmpfile)
  expect_identical(ncol(result), 5L)
})

# .validate_bed_ranges tests
test_that(".validate_bed_ranges accepts valid ranges", {
  bed <- data.frame(
    Chr = c("chr1", "chr1"),
    RegionStart = c(0L, 100L),
    RegionEnd = c(100L, 200L),
    stringsAsFactors = FALSE
  )
  expect_invisible(.validate_bed_ranges(bed))
})

test_that(".validate_bed_ranges rejects invalid ranges", {
  bed <- data.frame(
    Chr = c("chr1", "chr1"),
    RegionStart = c(0L, 100L),
    RegionEnd = c(100L, 100L),
    stringsAsFactors = FALSE
  )
  expect_error(
    .validate_bed_ranges(bed),
    regexp = "ERROR 4"
  )

  bed <- data.frame(
    Chr = "chr1",
    RegionStart = 200L,
    RegionEnd = 100L,
    stringsAsFactors = FALSE
  )
  expect_error(
    .validate_bed_ranges(bed),
    regexp = "Invalid BED"
  )
})

# .check_tab_file tests
test_that(".check_tab_file accepts valid files", {
  tmpfile <- tempfile()
  writeLines(
    c("chr1\t0\t100", "chr2\t100\t200"),
    tmpfile
  )
  on.exit(unlink(tmpfile))

  expect_invisible(
    .check_tab_file(tmpfile, min_fields = 3L, label = "BED", error_num = 1L)
  )
})

test_that(".check_tab_file skips comments", {
  tmpfile <- tempfile()
  writeLines(
    c("# Comment", "chr1\t0\t100"),
    tmpfile
  )
  on.exit(unlink(tmpfile))

  expect_invisible(
    .check_tab_file(tmpfile, min_fields = 3L, label = "BED", error_num = 1L)
  )
})

test_that(".check_tab_file validates field counts", {
  tmpfile <- tempfile()
  writeLines("chr1\t0", tmpfile)
  on.exit(unlink(tmpfile))

  expect_error(
    .check_tab_file(tmpfile, min_fields = 3L, label = "BED", error_num = 3L),
    regexp = "ERROR 3"
  )

  tmpfile2 <- tempfile()
  writeLines("chr1\t0\t100\textra\tfields", tmpfile2)
  on.exit(unlink(tmpfile2), add = TRUE)

  expect_error(
    .check_tab_file(
      tmpfile2,
      min_fields = 3L,
      max_fields = 3L,
      label = "BED",
      error_num = 3L
    ),
    regexp = "ERROR 3"
  )
})

# validate_inputs integration tests
test_that("validate_inputs returns valid structure", {
  gtf_tmp <- tempfile(fileext = ".gtf")
  bed_tmp <- tempfile(fileext = ".bed")

  writeLines(
    'chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1";',
    gtf_tmp
  )
  writeLines("chr1\t0\t100", bed_tmp)

  on.exit({
    unlink(gtf_tmp)
    unlink(bed_tmp)
  })

  result <- validate_inputs(gtf_tmp, bed_tmp)
  expect_type(result, "list")
  expect_named(result, c("gtf", "bed"))
  expect_s3_class(result$gtf, "data.frame")
  expect_s3_class(result$bed, "data.frame")
})

test_that("validate_inputs validates all inputs", {
  gtf_tmp <- tempfile(fileext = ".gtf")
  bed_tmp <- tempfile(fileext = ".bed")

  writeLines("chr1\tensembl\tgene", gtf_tmp) # Too few columns
  writeLines("chr1\t0\t100", bed_tmp)

  on.exit({
    unlink(gtf_tmp)
    unlink(bed_tmp)
  })

  expect_error(
    validate_inputs(gtf_tmp, bed_tmp),
    regexp = "ERROR 2"
  )

  # GTF valid, BED invalid
  writeLines(
    'chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1";',
    gtf_tmp
  )
  writeLines("chr1\t100\t100", bed_tmp) # Invalid range

  expect_error(
    validate_inputs(gtf_tmp, bed_tmp),
    regexp = "ERROR 4"
  )
})
