test_that("process_gtf works on example files and preserves order", {
  gtf_file <- system.file(
    "extdata",
    "AlgorithmToy.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  reg_file <- system.file(
    "extdata",
    "AlgorithmToy.bed",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  out <- process_gtf(reg_file, gtf_file, genome_name = "test", genome_version = "v1")

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
  expect_true(all(grepl(":", out$Chr, fixed = TRUE)))
  expect_true(all(grepl("-", out$Chr, fixed = TRUE)))

  # Start/end must be positive after shifting
  expect_true(all(out$start >= 1L))
  expect_true(all(out$end >= out$start))

  # Output row order must match input GTF order
  gtf_raw <- read.table(
    gtf_file,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE
  )

  # attribute is column 9
  expect_identical(nrow(out), nrow(gtf_raw))
  expect_identical(out$attribute, gtf_raw[[9]])
})


test_that("process_gtf aborts when split regions overlap (ambiguous)", {
  gtf_file <- system.file(
    "extdata",
    "AlgorithmToy.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  # Load BED using base R
  reg <- read.table(
    system.file(
      "extdata",
      "AlgorithmToy.bed",
      package = "scShardSplitRef",
      mustWork = TRUE
    ),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(reg)[1:3] <- c("Chr", "RegionStart", "RegionEnd")

  # Force overlap on chr1
  reg$RegionStart[1] <- reg$RegionStart[2] - 5
  reg$RegionEnd[1] <- reg$RegionEnd[2] + 5

  reg_file <- tempfile(fileext = ".bed")
  write.table(
    reg,
    reg_file,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  expect_error(
    process_gtf(reg_file, gtf_file,  genome_name = "test", genome_version = "v1"),
    "Ambiguous split-region matches"
  )
})


test_that("process_gtf aborts when some features match no region", {
  gtf_file <- system.file(
    "extdata",
    "AlgorithmToy.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  # Load BED using base R
  reg <- read.table(
    system.file(
      "extdata",
      "AlgorithmToy.bed",
      package = "scShardSplitRef",
      mustWork = TRUE
    ),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(reg)[1:3] <- c("Chr", "RegionStart", "RegionEnd")

  # Remove all chr2 regions
  reg <- reg[reg$Chr != "chr2", , drop = FALSE]

  reg_file <- tempfile(fileext = ".bed")
  write.table(
    reg,
    reg_file,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )

  expect_error(
    process_gtf(reg_file, gtf_file, genome_name = "test", genome_version = "v1"),
    "do not fall within any split-region interval"
  )
})


test_that(".prepare_split_regions constructs NEW correctly", {
  regions <- data.frame(
    Chr = "chrX",
    RegionStart = 100,
    RegionEnd = 200,
    stringsAsFactors = FALSE
  )

  out <- scShardSplitRef:::`.prepare_split_regions`(regions)

  expect_identical(out$NEW, "chrX:100-200")
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
    c(FALSE, TRUE, TRUE, TRUE)
  )
})
