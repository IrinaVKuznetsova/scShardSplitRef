test_that("process_gtf works on example files and preserves order, writes a correctly ordered GTF to disk", {
  gtf_file <- system.file(
    "extdata", "AlgorithmToy.gtf",
    package = "scShardSplitRef", mustWork = TRUE
  )
  reg_file <- system.file(
    "outdata", "DetermineSplitReg_AlgorithmToy.bed",
    package = "scShardSplitRef", mustWork = TRUE
  )
  
  genome_name    <- "synthetic"
  genome_version <- "v1"
  keep_attributes <- NULL   ## (FALSE was wrong; NULL - "keep all attributes")
  out_path       <- tempdir()
  
  res <- process_gtf(reg_file, gtf_file, genome_name, genome_version,
                     keep_attributes=NULL, out_path)
  
  # Write to disk and return NULL (in the initial code this part was a data.frame, thus checks were failing now)
  expect_null(res)
  
  # write output to disk, and read back in the next line so all checks can be done
  written <- file.path(out_path,
                       "B1_FINAL_MODIFIED_GTF_synthetic_v1.gtf")
  expect_true(file.exists(written))
  
  # Read the just written file back and run checks
  out <- read.table(written, sep = "\t", header = FALSE, quote = "",
                    stringsAsFactors = FALSE)
  colnames(out) <- c("Chr", "source", "feature", "start", "end",
                     "score", "strand", "frame", "attribute")
  
  gtf_raw <- read.table(gtf_file, sep = "\t", header = FALSE, quote = "",
                        comment.char = "#", stringsAsFactors = FALSE)
  
  # Same number of rows, in the same order
  expect_identical(nrow(out), nrow(gtf_raw))
  
  # Chr column has split-region e.g. "chr1:0-410"
  expect_true(all(grepl(":", out$Chr, fixed = TRUE)))
  expect_true(all(grepl("-", out$Chr, fixed = TRUE)))
  
  # Coordinates are valid after shifting
  expect_true(all(out$start >= 1L))
  expect_true(all(out$end >= out$start))
  
  # Attributes are unchanged (keep_attributes = NULL)
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
      "outdata",
      "DetermineSplitReg_AlgorithmToy.bed",
      package = "scShardSplitRef",
      mustWork = TRUE
    ),
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE
  )

  colnames(reg)[1:3] <- c("Chr", "RegionStart", "RegionEnd")

  # Force overlap on chr1 - hard coded
  #reg$RegionStart[1] <- reg$RegionStart[2] - 5
  #reg$RegionEnd[1] <- reg$RegionEnd[2] + 5
  chr1_rows <- which(reg$Chr == "chr1")   # removing hard coded lines, changed to chr1 
  reg$RegionStart[chr1_rows[2]] <- reg$RegionStart[chr1_rows[1]]  # [2/1]chr regions Start and End

  reg_file <- tempfile(fileext = ".bed")
  write.table(
    reg,
    reg_file,
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
  
  # Define param
  genome_name     <- "synthetic"
  genome_version  <- "v1"
  keep_attributes <- NULL
  out_path        <- tempdir()

  expect_error(
    process_gtf(reg_file, gtf_file, genome_name, genome_version, keep_attributes=NULL, out_path),
    "ambiguous split-region matches"
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
      "outdata",
      "DetermineSplitReg_AlgorithmToy.bed",
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

  # Define param
  genome_name     <- "synthetic"
  genome_version  <- "v1"
  keep_attributes <- NULL
  out_path        <- tempdir()
  
  expect_error(
    process_gtf(reg_file, gtf_file, genome_name, genome_version, keep_attributes=NULL, out_path),
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
