# Tests ensuring genomic coordinates remain true integers throughout the
# pipeline, and are never written to disk in scientific notation
# (e.g. "1e+08" instead of "100000000"). Covers the two defensive helpers
# (.coerce_integer_coord in validate_inputs.R, .format_coord in
# process_gtf.R) directly, plus end-to-end checks on the actual files
# written by determine_split_regions() and process_gtf().

# .coerce_integer_coord ---------------------------------------------------------

test_that(".coerce_integer_coord returns integer type for whole-number doubles", {
  result <- .coerce_integer_coord(c(0, 100000000), "start")
  expect_type(result, "integer")
  expect_identical(result, c(0L, 100000000L))
})

test_that(".coerce_integer_coord passes already-integer input through unchanged", {
  x <- c(1L, 2L, 3L)
  expect_identical(.coerce_integer_coord(x, "start"), x)
})

test_that(".coerce_integer_coord coerces numeric-looking character input", {
  result <- .coerce_integer_coord(c("100", "200000000"), "start")
  expect_type(result, "integer")
  expect_identical(result, c(100L, 200000000L))
})

test_that(".coerce_integer_coord errors on fractional (non-whole) values", {
  expect_error(
    .coerce_integer_coord(c(1, 2.5, 3), "start"),
    regexp = "whole-number coordinates"
  )
})

test_that(".coerce_integer_coord errors on values outside the 32-bit integer range", {
  too_big <- .Machine$integer.max + 1
  expect_error(
    .coerce_integer_coord(too_big, "start"),
    regexp = "32-bit integer range"
  )
})

test_that(".coerce_integer_coord does not silently produce NA on overflow", {
  # Guard against the exact footgun this helper exists to prevent: plain
  # as.integer() would return NA + a warning here instead of erroring.
  too_big <- .Machine$integer.max * 2
  expect_error(.coerce_integer_coord(too_big, "start"))
})

# .format_coord -------------------------------------------------------------

test_that(".format_coord renders large round doubles without scientific notation", {
  # 1e8 is the classic R footgun: as.character(1e8) == "1e+08"
  expect_identical(.format_coord(1e8), "100000000")
  expect_identical(.format_coord(1e9), "1000000000")
  expect_identical(.format_coord(100000000L), "100000000")
})

test_that(".format_coord handles a vector mixing round and non-round values", {
  result <- .format_coord(c(1e8, 100000001, 5))
  expect_identical(result, c("100000000", "100000001", "5"))
})

test_that(".format_coord passes NA through as NA, not the string 'NA'", {
  result <- .format_coord(c(1e8, NA, 5))
  expect_identical(result, c("100000000", NA_character_, "5"))
})

test_that(".format_coord handles zero and negative coordinates", {
  expect_identical(.format_coord(0), "0")
  expect_identical(.format_coord(-5), "-5")
})

test_that(".format_coord errors on non-numeric input", {
  expect_error(.format_coord("100000000"), regexp = "must be numeric")
  expect_error(.format_coord(TRUE), regexp = "must be numeric")
})

# .build_boundaries / .subdivide_oversized stay integer at large scale --------

test_that(".build_boundaries returns integer boundaries for near-chromosome-scale input", {
  # ~2 Gb interval, comfortably inside 32-bit range but large enough that
  # naive double arithmetic without a final as.integer() would round-trip
  # to scientific notation on write.
  result <- .build_boundaries(0L, 2000000000L, 536870912L)
  expect_type(result, "integer")
  expect_false(anyNA(result))
  expect_true(all(result > 0L & result < 2000000000L))
})

test_that(".subdivide_oversized returns integer boundaries for large intervals", {
  result <- .subdivide_oversized(c(0L, 1000000000L), limit = 536870912L)
  expect_type(result, "integer")
  expect_false(anyNA(result))
})

# .process_single_region correctly avoids a real gene ---------------------
#
# NOTE: the existing tests for .process_single_region() (in
# test-determine_split_regions.R / test-split_regions.R) all pass an empty
# data.frame as the `merged_gene_ranges` argument, but the function actually
# expects a *named list keyed by chromosome* (as produced by
# .merge_all_gene_ranges()). Because those tests use zero genes, indexing
# the data.frame with `df[["chr1"]]` happens to silently return NULL and
# fall back to "no genes" -- so the tests pass, but they never exercise gene
# avoidance at all. This test uses the correctly-shaped input.

test_that(".process_single_region shifts a boundary that would fall inside a gene", {
  region <- data.frame(
    chr = "chr1",
    start = 0L,
    end = 300L,
    stringsAsFactors = FALSE
  )

  # A gene straddling the natural split point at width/2 = 150.
  merged_gene_ranges <- list(
    chr1 = data.frame(start = 140L, end = 160L)
  )

  result <- .process_single_region(
    region,
    merged_gene_ranges,
    limit = 150L,
    shift_by = 1L,
    clearance = 1L
  )

  boundaries <- result$end[-nrow(result)]
  expect_true(all(boundaries <= 139L | boundaries >= 161L))
})

# End-to-end: determine_split_regions() never writes scientific notation ------

test_that("determine_split_regions() writes large coordinates without scientific notation", {
  bed <- withr::local_tempfile(fileext = ".bed")
  gtf <- withr::local_tempfile(fileext = ".gtf")
  out <- withr::local_tempfile(fileext = ".bed")

  # A single ~1 Gb region, forcing the size-limit splitting logic to
  # generate several large, round-ish boundary values.
  writeLines("chr1\t0\t1000000000", bed)
  writeLines(
    'chr1\tsrc\tgene\t500000000\t500000100\t.\t+\t.\tgene_id "G1";',
    gtf
  )

  determine_split_regions(
    bed = bed,
    gtf = gtf,
    output_bed = out,
    limit = 2L^29L,
    shift_by = 1L,
    clearance = 1L
  )

  written_lines <- readLines(out)
  expect_gt(length(written_lines), 1L)
  expect_false(any(grepl(
    "[0-9]e[+-]?[0-9]",
    written_lines,
    ignore.case = TRUE
  )))

  parsed <- read.table(out, sep = "\t", stringsAsFactors = FALSE)
  colnames(parsed) <- c("chr", "start", "end")
  expect_type(parsed$start, "integer")
  expect_type(parsed$end, "integer")

  # Intervals tile the original region exactly, with no gaps or overlaps,
  # and none exceed the requested limit.
  ord <- order(parsed$start)
  parsed <- parsed[ord, ]
  expect_identical(parsed$start[1], 0L)
  expect_identical(parsed$end[nrow(parsed)], 1000000000L)
  expect_identical(parsed$start[-1], parsed$end[-nrow(parsed)])
  expect_true(all(parsed$end - parsed$start <= 2L^29L))
})

# End-to-end: process_gtf() never writes scientific notation ------------------

test_that("process_gtf() writes large shifted coordinates without scientific notation", {
  reg_file <- withr::local_tempfile(fileext = ".bed")
  gtf_file <- withr::local_tempfile(fileext = ".gtf")
  out_path <- withr::local_tempdir()

  # RegionStart = 0 and a feature chosen so the shifted coordinate
  # (start - RegionStart + 1) lands exactly on 1e8/1e8+100 -- the classic
  # round numbers that silently render in scientific notation without
  # .format_coord().
  writeLines("chr1\t0\t250000000", reg_file)
  writeLines(
    'chr1\tsrc\tgene\t99999999\t100000099\t.\t+\t.\tgene_id "G1";',
    gtf_file
  )

  process_gtf(
    split_regions_bed = reg_file,
    gtf = gtf_file,
    genome_name = "bigcoord",
    genome_version = "v1",
    out_path = out_path
  )

  written_file <- file.path(out_path, "B1_FINAL_MODIFIED_GTF_bigcoord_v1.gtf")
  expect_true(file.exists(written_file))

  raw_lines <- readLines(written_file)
  expect_false(any(grepl("[0-9]e[+-]?[0-9]", raw_lines, ignore.case = TRUE)))

  out <- read.table(written_file, sep = "\t", stringsAsFactors = FALSE)
  colnames(out) <- c(
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

  expect_identical(out$start, 100000000L)
  expect_identical(out$end, 100000100L)
})
