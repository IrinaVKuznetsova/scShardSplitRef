test_that("process_gtf works on example files and preserves order", {
  gtf_file <- system.file(
    "extdata",
    "A3_toy_all_scenarios_2chr.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  reg_file <- system.file(
    "extdata",
    "IN0_toy_centromeres_for_gtf.bed",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  out <- process_gtf(gtf_file, reg_file)

  # Basic structure checks
  expect_s3_class(out, "data.frame")
  expect_true(all(
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
    ) %in%
      colnames(out)
  ))

  # Chr column must contain split-region names
  expect_true(all(grepl(":", out$Chr)))
  expect_true(all(grepl("-", out$Chr)))

  # Start/end must be positive after shifting
  expect_true(all(out$start >= 1))
  expect_true(all(out$end >= out$start))

  # Output row order must match input GTF order
  gtf_raw <- readr::read_tsv(gtf_file, show_col_types = FALSE)
  expect_equal(nrow(out), nrow(gtf_raw))
  expect_equal(out$attribute, gtf_raw$attribute)
})


test_that("process_gtf aborts when split regions overlap (ambiguous)", {
  gtf_file <- system.file(
    "extdata",
    "A3_toy_all_scenarios_2chr.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  # Load example BED and artificially create overlapping regions
  reg <- readr::read_tsv(
    system.file(
      "extdata",
      "IN0_toy_centromeres_for_gtf.bed",
      package = "scShardSplitRef",
      mustWork = TRUE
    ),
    show_col_types = FALSE
  )

  # Force overlap on chr1
  reg$RegionStart[1] <- reg$RegionStart[2] - 5
  reg$RegionEnd[1] <- reg$RegionEnd[2] + 5

  reg_file <- tempfile(fileext = ".bed")
  readr::write_tsv(reg, reg_file)

  expect_error(
    process_gtf(gtf_file, reg_file),
    "Ambiguous split-region matches"
  )
})


test_that("process_gtf aborts when some features match no region", {
  gtf_file <- system.file(
    "extdata",
    "A3_toy_all_scenarios_2chr.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  # Load example BED and remove all chr2 regions
  reg <- readr::read_tsv(
    system.file(
      "extdata",
      "IN0_toy_centromeres_for_gtf.bed",
      package = "scShardSplitRef",
      mustWork = TRUE
    ),
    show_col_types = FALSE
  )

  reg <- reg[reg$Chr != "chr2", ]

  reg_file <- tempfile(fileext = ".bed")
  readr::write_tsv(reg, reg_file)

  expect_error(
    process_gtf(gtf_file, reg_file),
    "do not fit within any split-region interval"
  )
})


test_that(".prepare_split_regions constructs NEW correctly", {
  regions <- data.frame(
    Chr = "chrX",
    RegionStart = 100,
    RegionEnd = 200
  )

  out <- scShardSplitRef:::`.prepare_split_regions`(regions)

  expect_equal(out$NEW, "chrX:100-200")
})

test_that(".rows_within_region applies matching rules correctly", {
  df <- data.frame(
    start = c(5, 10, 10, 11),
    end = c(8, 12, 10, 15),
    RegionStart = 10,
    RegionEnd = 20
  )

  expect_identical(
    scShardSplitRef:::`.rows_within_region`(df),
    c(FALSE, TRUE, TRUE, FALSE)
  )
})
