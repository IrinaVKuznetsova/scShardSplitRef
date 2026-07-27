# .find_invalid_field_counts --------------------------------------------------

test_that(".find_invalid_field_counts returns empty for valid counts", {
  expect_identical(
    .find_invalid_field_counts(c(3L, 3L, 4L), min_fields = 3L),
    integer(0)
  )
})

test_that(".find_invalid_field_counts catches below min", {
  expect_identical(
    .find_invalid_field_counts(c(3L, 2L, 3L), min_fields = 3L),
    2L
  )
})

test_that(".find_invalid_field_counts catches above max", {
  expect_identical(
    .find_invalid_field_counts(c(3L, 5L, 3L), min_fields = 3L, max_fields = 4L),
    2L
  )
})

test_that(".find_invalid_field_counts catches both below and above bounds", {
  result <- .find_invalid_field_counts(
    c(2L, 3L, 6L),
    min_fields = 3L,
    max_fields = 4L
  )
  expect_identical(result, c(1L, 3L))
})

test_that(".find_invalid_field_counts deduplicates and sorts", {
  result <- .find_invalid_field_counts(
    c(5L, 1L, 5L),
    min_fields = 3L,
    max_fields = 4L
  )
  expect_identical(result, c(1L, 2L, 3L)) # 1 & 3 exceed max, 2 below min
})

# .format_field_expectation ---------------------------------------------------

test_that(".format_field_expectation with no max", {
  expect_identical(
    .format_field_expectation(9L, NA_integer_),
    "at least 9 TAB-separated fields"
  )
})

test_that(".format_field_expectation with min and max", {
  expect_identical(
    .format_field_expectation(3L, 9L),
    "3-9 TAB-separated fields"
  )
})

# .check_files_exist ----------------------------------------------------------

test_that(".check_files_exist passes when both files exist", {
  bed <- withr::local_tempfile(fileext = ".bed")
  gtf <- withr::local_tempfile(fileext = ".gtf")
  file.create(bed)
  file.create(gtf)

  expect_invisible(.check_files_exist(bed, gtf))
})

test_that(".check_files_exist errors when BED missing", {
  gtf <- withr::local_tempfile(fileext = ".gtf")
  file.create(gtf)

  expect_error(
    .check_files_exist("nonexistent.bed", gtf),
    "BED file does not exist"
  )
})

test_that(".check_files_exist errors when GTF missing", {
  bed <- withr::local_tempfile(fileext = ".bed")
  file.create(bed)

  expect_error(
    .check_files_exist(bed, "nonexistent.gtf"),
    "GTF file does not exist"
  )
})

# .validate_bed_ranges --------------------------------------------------------

test_that(".validate_bed_ranges passes for valid ranges", {
  bed <- data.frame(RegionStart = c(100L, 300L), RegionEnd = c(200L, 400L))
  expect_invisible(.validate_bed_ranges(bed))
})

test_that(".validate_bed_ranges errors when start == end", {
  bed <- data.frame(RegionStart = c(100L, 200L), RegionEnd = c(200L, 200L))
  expect_error(.validate_bed_ranges(bed), "Invalid BED coordinates")
})

test_that(".validate_bed_ranges errors when start > end", {
  bed <- data.frame(RegionStart = c(500L), RegionEnd = c(100L))
  expect_error(.validate_bed_ranges(bed), "Invalid BED coordinates")
})

test_that(".validate_bed_ranges error includes problematic row number", {
  bed <- data.frame(
    RegionStart = c(100L, 500L, 300L),
    RegionEnd = c(200L, 100L, 400L)
  )
  expect_error(.validate_bed_ranges(bed), "2")
})

# validate_inputs (integration) -----------------------------------------------

test_that("validate_inputs returns list with bed and gtf", {
  bed_tmp <- withr::local_tempfile(fileext = ".bed")
  gtf_tmp <- withr::local_tempfile(fileext = ".gtf")

  writeLines(
    c(
      "chr1\t100\t200",
      "chr1\t300\t400"
    ),
    bed_tmp
  )

  writeLines(
    c(
      "chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute",
      "chr1\tsource\texon\t1\t50\t.\t+\t.\tattribute"
    ),
    gtf_tmp
  )

  result <- .validate_inputs(bed_tmp, gtf_tmp)

  expect_named(result, c("bed", "gtf"))
  expect_s3_class(result$bed, "data.frame")
  expect_s3_class(result$gtf, "data.frame")
})

test_that("validate_inputs assigns correct BED column names", {
  bed_tmp <- withr::local_tempfile(fileext = ".bed")
  gtf_tmp <- withr::local_tempfile(fileext = ".gtf")

  writeLines("chr1\t100\t200", bed_tmp)
  writeLines("chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute", gtf_tmp)

  result <- .validate_inputs(bed_tmp, gtf_tmp)

  expect_named(result$bed[1:3], c("Chr", "RegionStart", "RegionEnd"))
})

test_that("validate_inputs assigns correct GTF column names", {
  bed_tmp <- withr::local_tempfile(fileext = ".bed")
  gtf_tmp <- withr::local_tempfile(fileext = ".gtf")

  writeLines("chr1\t100\t200", bed_tmp)
  writeLines("chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute", gtf_tmp)

  result <- .validate_inputs(bed_tmp, gtf_tmp)

  expect_named(
    result$gtf,
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

test_that("validate_inputs errors on invalid BED ranges", {
  bed_tmp <- withr::local_tempfile(fileext = ".bed")
  gtf_tmp <- withr::local_tempfile(fileext = ".gtf")

  writeLines("chr1\t500\t100", bed_tmp) # start > end
  writeLines("chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute", gtf_tmp)

  expect_error(.validate_inputs(bed_tmp, gtf_tmp), "Invalid BED coordinates")
})
