test_that(".check_tab_file passes for valid BED file", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr1\t100\t200",
      "chr1\t300\t400"
    ),
    tmp
  )

  expect_invisible(.check_tab_file(tmp, min_fields = 3L, label = "BED"))
})

test_that(".check_tab_file passes for valid GTF file", {
  tmp <- withr::local_tempfile(fileext = ".gtf")
  writeLines(
    "chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute",
    tmp
  )

  expect_invisible(.check_tab_file(tmp, min_fields = 9L, label = "GTF"))
})

test_that(".check_tab_file skips comment and blank lines", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(
    c(
      "# this is a comment",
      "",
      "   ",
      "chr1\t100\t200"
    ),
    tmp
  )

  expect_invisible(.check_tab_file(tmp, min_fields = 3L, label = "BED"))
})

test_that(".check_tab_file returns invisible TRUE for all-comment file", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(c("# comment only", "# another comment"), tmp)

  expect_invisible(.check_tab_file(tmp, min_fields = 3L, label = "BED"))
  expect_true(.check_tab_file(tmp, min_fields = 3L, label = "BED"))
})

test_that(".check_tab_file errors on wrong file extension", {
  tmp <- withr::local_tempfile(fileext = ".txt")
  writeLines("chr1\t100\t200", tmp)

  expect_error(
    .check_tab_file(tmp, min_fields = 3L, label = "BED"),
    "File extension does not match"
  )
})

test_that(".check_tab_file errors on non-tab-separated file", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr1,100,200",
      "chr1,300,400"
    ),
    tmp
  )

  expect_error(
    .check_tab_file(tmp, min_fields = 3L, label = "BED"),
    regexp = "tab-separated",
    ignore.case = TRUE
  )
})

test_that(".check_tab_file errors when too few fields", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr1\t100", # only 2 fields, need 3
      "chr1\t300"
    ),
    tmp
  )

  expect_error(
    .check_tab_file(tmp, min_fields = 3L, label = "BED"),
    "Invalid BED file format"
  )
})

test_that(".check_tab_file errors when too many fields and max_fields set", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr1\t100\t200\textra\tfield"
    ),
    tmp
  )

  expect_error(
    .check_tab_file(tmp, min_fields = 3L, max_fields = 3L, label = "BED"),
    "Invalid BED file format"
  )
})

test_that(".check_tab_file error message includes problematic line number", {
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(
    c(
      "chr1\t100\t200",
      "chr1\t300", # line 2 is bad
      "chr1\t400\t500"
    ),
    tmp
  )

  expect_error(
    .check_tab_file(tmp, min_fields = 3L, label = "BED"),
    "2"
  )
})

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

  result <- validate_inputs(bed_tmp, gtf_tmp)

  expect_named(result, c("bed", "gtf"))
  expect_s3_class(result$bed, "data.frame")
  expect_s3_class(result$gtf, "data.frame")
})

test_that("validate_inputs assigns correct BED column names", {
  bed_tmp <- withr::local_tempfile(fileext = ".bed")
  gtf_tmp <- withr::local_tempfile(fileext = ".gtf")

  writeLines("chr1\t100\t200", bed_tmp)
  writeLines("chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute", gtf_tmp)

  result <- validate_inputs(bed_tmp, gtf_tmp)

  expect_named(result$bed[1:3], c("Chr", "RegionStart", "RegionEnd"))
})

test_that("validate_inputs assigns correct GTF column names", {
  bed_tmp <- withr::local_tempfile(fileext = ".bed")
  gtf_tmp <- withr::local_tempfile(fileext = ".gtf")

  writeLines("chr1\t100\t200", bed_tmp)
  writeLines("chr1\tsource\tgene\t1\t100\t.\t+\t.\tattribute", gtf_tmp)

  result <- validate_inputs(bed_tmp, gtf_tmp)

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

  expect_error(validate_inputs(bed_tmp, gtf_tmp), "Invalid BED coordinates")
})
