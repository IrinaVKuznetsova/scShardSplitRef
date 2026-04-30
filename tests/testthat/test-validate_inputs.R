# .check_files_exist tests
test_that(".check_files_exist validates file existence", {
  expect_invisible(.check_files_exist(
    system.file(
      "extdata",
      "A3_toy_all_scenarios_2chr.gtf",
      package = "scShardSplitRef"
    ),
    system.file(
      "extdata",
      "IN0_toy_centromeres_for_gtf.bed",
      package = "scShardSplitRef"
    )
  ))
})

test_that(".check_files_exist reports missing GTF", {
  expect_error(
    .check_files_exist("/nonexistent/file.gtf", tempfile()),
    regexp = "GTF file does not exist"
  )
})

test_that(".check_files_exist reports missing BED", {
  gtf_path <- system.file(
    "extdata",
    "A3_toy_all_scenarios_2chr.gtf",
    package = "scShardSplitRef"
  )
  expect_error(
    .check_files_exist(gtf_path, "/nonexistent/file.bed"),
    regexp = "BED file does not exist"
  )
})

# .find_invalid_field_counts tests
test_that(".find_invalid_field_counts identifies too few fields", {
  result <- .find_invalid_field_counts(c(2, 3, 2, 5), 3L, NA_integer_)
  expect_identical(result, c(1L, 3L))
})

test_that(".find_invalid_field_counts identifies too many fields", {
  result <- .find_invalid_field_counts(c(3, 6, 3, 8), 3L, 5L)
  expect_identical(result, c(2L, 4L))
})

test_that(".find_invalid_field_counts returns empty for valid counts", {
  result <- .find_invalid_field_counts(c(3, 4, 5), 3L, 5L)
  expect_identical(result, integer(0))
})

# .format_field_expectation tests
test_that(".format_field_expectation formats range correctly", {
  result <- .format_field_expectation(3L, 9L)
  expect_identical(result, "3-9 TAB-separated fields")
})

test_that(".format_field_expectation formats minimum only", {
  result <- .format_field_expectation(9L, NA_integer_)
  expect_identical(result, "at least 9 TAB-separated fields")
})

# .read_and_format_gtf tests
test_that(".read_and_format_gtf reads and formats", {
  tmpfile <- tempfile(fileext = ".gtf")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "chr1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\";",
      "chr1\tensembl\texon\t10\t50\t.\t+\t.\tgene_id \"G1\"; exon_id \"E1\";"
    ),
    tmpfile
  )

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
  expect_identical(nrow(result), 2L)
})

test_that(".read_and_format_gtf skips comments", {
  tmpfile <- tempfile(fileext = ".gtf")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "# This is a comment",
      "chr1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\";",
      "# Another comment",
      "chr1\tensembl\texon\t10\t50\t.\t+\t.\tgene_id \"G1\"; exon_id \"E1\";"
    ),
    tmpfile
  )

  result <- .read_and_format_gtf(tmpfile)
  expect_identical(nrow(result), 2L)
})

# .read_and_format_bed tests
test_that(".read_and_format_bed reads and formats", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "chr1\t0\t1000",
      "chr2\t500\t2000"
    ),
    tmpfile
  )

  result <- .read_and_format_bed(tmpfile)

  expect_identical(colnames(result)[1:3], c("Chr", "RegionStart", "RegionEnd"))
  expect_identical(nrow(result), 2L)
  expect_identical(result$RegionStart, c(0L, 500L)) # Use integers
})

test_that(".read_and_format_bed preserves extra columns", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "chr1\t0\t1000\textra1\textra2",
      "chr2\t500\t2000\textra3\textra4"
    ),
    tmpfile
  )

  result <- .read_and_format_bed(tmpfile)

  expect_identical(ncol(result), 5L)
  expect_identical(colnames(result)[4:5], c("V4", "V5"))
})

# .validate_bed_ranges tests
test_that(".validate_bed_ranges accepts valid ranges", {
  bed <- data.frame(
    Chr = c("chr1", "chr2"),
    RegionStart = c(0L, 100L),
    RegionEnd = c(1000L, 2000L),
    stringsAsFactors = FALSE
  )

  expect_invisible(.validate_bed_ranges(bed))
})

test_that(".validate_bed_ranges rejects invalid ranges", {
  bed <- data.frame(
    Chr = c("chr1", "chr2"),
    RegionStart = c(0L, 2000L),
    RegionEnd = c(1000L, 100L),
    stringsAsFactors = FALSE
  )

  expect_error(
    .validate_bed_ranges(bed),
    regexp = "Invalid BED coordinates"
  )
})

test_that(".validate_bed_ranges rejects equal start and end", {
  bed <- data.frame(
    Chr = "chr1",
    RegionStart = 100L,
    RegionEnd = 100L,
    stringsAsFactors = FALSE
  )

  expect_error(
    .validate_bed_ranges(bed),
    regexp = "Invalid BED coordinates"
  )
})

# .check_tab_file tests
test_that(".check_tab_file accepts valid files", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "chr1\t0\t1000",
      "chr2\t500\t2000"
    ),
    tmpfile
  )

  expect_invisible(.check_tab_file(tmpfile, min_fields = 3L, label = "BED"))
})

test_that(".check_tab_file skips comments", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "# Comment line",
      "chr1\t0\t1000",
      "# Another comment",
      "chr2\t500\t2000"
    ),
    tmpfile
  )

  expect_invisible(.check_tab_file(tmpfile, min_fields = 3L, label = "BED"))
})

test_that(".check_tab_file validates field counts", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  writeLines(
    c(
      "chr1\t0\t1000",
      "chr2\t500",
      "chr3\t100\t200\t300"
    ),
    tmpfile
  )

  expect_error(
    .check_tab_file(tmpfile, min_fields = 3L, label = "BED"),
    regexp = "Invalid BED file format"
  )
})

# validate_inputs integration tests
test_that("validate_inputs validates all inputs", {
  skip_if_not_installed("cli")

  gtf_tmp <- tempfile(fileext = ".gtf")
  bed_tmp <- tempfile(fileext = ".bed")
  on.exit({
    unlink(gtf_tmp)
    unlink(bed_tmp)
  })

  # Create invalid GTF (too few columns)
  writeLines("chr1 ensembl gene", gtf_tmp)
  writeLines("chr1\t0\t1000", bed_tmp)

  expect_error(
    validate_inputs(gtf_tmp, bed_tmp),
    regexp = "Invalid GTF file format"
  )
})

test_that("validate_inputs returns correct structure", {
  skip_if_not_installed("cli")

  gtf_tmp <- tempfile(fileext = ".gtf")
  bed_tmp <- tempfile(fileext = ".bed")
  on.exit({
    unlink(gtf_tmp)
    unlink(bed_tmp)
  })

  writeLines(
    c(
      "chr1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\";",
      "chr1\tensembl\texon\t10\t50\t.\t+\t.\tgene_id \"G1\"; exon_id \"E1\";"
    ),
    gtf_tmp
  )
  writeLines("chr1\t0\t1000", bed_tmp)

  result <- validate_inputs(gtf_tmp, bed_tmp)

  expect_type(result, "list")
  expect_true("gtf" %in% names(result))
  expect_true("bed" %in% names(result))
  expect_identical(nrow(result$gtf), 2L)
  expect_identical(nrow(result$bed), 1L)
})

test_that("validate_inputs rejects invalid BED ranges", {
  skip_if_not_installed("cli")

  gtf_tmp <- tempfile(fileext = ".gtf")
  bed_tmp <- tempfile(fileext = ".bed")
  on.exit({
    unlink(gtf_tmp)
    unlink(bed_tmp)
  })

  writeLines(
    "chr1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id \"G1\";",
    gtf_tmp
  )
  writeLines("chr1\t1000\t500", bed_tmp)

  expect_error(
    validate_inputs(gtf_tmp, bed_tmp),
    regexp = "Invalid BED coordinates"
  )
})
