# Helper to create test gene ranges
.make_gene_ranges <- function(starts, ends) {
  data.frame(
    start = starts,
    end = ends,
    stringsAsFactors = FALSE
  )
}

# Helper to create test GTF data
.make_gtf <- function(chrs, starts, ends, ids) {
  data.frame(
    chr = chrs,
    start = starts,
    end = ends,
    unique_ID = ids,
    stringsAsFactors = FALSE
  )
}

# .build_boundaries tests
test_that(".build_boundaries respects size limit", {
  expect_identical(.build_boundaries(0L, 100L, 200L), integer(0L))
  result <- .build_boundaries(0L, 1000L, 250L)
  expect_type(result, "integer")
  expect_true(length(result) > 0L)
})

test_that(".build_boundaries creates valid chunks", {
  boundaries <- .build_boundaries(0L, 1000L, 250L)
  points <- c(0L, boundaries, 1000L)
  widths <- diff(points)
  expect_true(all(widths <= 250L))
  expect_true(all(widths > 0L))
})

test_that(".build_boundaries handles exact multiples", {
  boundaries <- .build_boundaries(0L, 400L, 100L)
  points <- c(0L, boundaries, 400L)
  expect_identical(sum(diff(points)), 400L)
})

# .validate_bed_coordinates tests
test_that(".validate_bed_coordinates accepts valid coords", {
  expect_invisible(.validate_bed_coordinates(0L, 100L, "chr1"))
  expect_invisible(.validate_bed_coordinates(1000L, 2000L, "chr1"))
})

test_that(".validate_bed_coordinates rejects invalid coords", {
  expect_error(
    .validate_bed_coordinates(-1L, 100L, "chr1"),
    regexp = ">= 0"
  )
  expect_error(
    .validate_bed_coordinates(100L, 100L, "chr1"),
    regexp = "must be < end"
  )
  expect_error(
    .validate_bed_coordinates(100L, 100L, "chrX"),
    regexp = "chrX"
  )
})

# .boundaries_to_df tests
test_that(".boundaries_to_df creates correct structure", {
  df <- .boundaries_to_df("chr1", 0L, 100L, c(30L, 60L))
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), 3L)
  expect_identical(colnames(df), c("chr", "start", "end"))
  expect_identical(df$start, c(0L, 30L, 60L))
  expect_identical(df$end, c(30L, 60L, 100L))
})

test_that(".boundaries_to_df handles edge cases", {
  # Empty boundaries
  df <- .boundaries_to_df("chr1", 0L, 50L, integer(0L))
  expect_identical(nrow(df), 1L)
  expect_identical(df$start, 0L)
  expect_identical(df$end, 50L)

  # Single boundary
  df <- .boundaries_to_df("chrX", 0L, 50L, 25L)
  expect_identical(df$chr, c("chrX", "chrX"))
})

test_that(".boundaries_to_df ensures non-overlapping intervals", {
  df <- .boundaries_to_df("chr1", 0L, 1000L, c(250L, 500L, 750L))
  for (i in seq_len(nrow(df) - 1)) {
    expect_identical(df$end[i], df$start[i + 1])
  }
})

# .split_region_core tests
test_that(".split_region_core respects size limit", {
  result <- .split_region_core("chr1", 0L, 1000L, 250L)
  expect_s3_class(result, "data.frame")
  expect_true(all(result$end - result$start <= 250L))

  result <- .split_region_core("chr1", 0L, 100L, 500L)
  expect_identical(nrow(result), 1L)
  expect_identical(result$start, 0L)
})

test_that(".split_region_core preserves chromosome", {
  result <- .split_region_core("chrY", 0L, 100L, 50L)
  expect_identical(unique(result$chr), "chrY")
})

# .merge_overlapping_gene_ranges tests
test_that(".merge_overlapping_gene_ranges merges overlaps", {
  genes <- .make_gene_ranges(c(10L, 15L), c(30L, 25L))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 1L)
  expect_identical(merged$start, 10L)
  expect_identical(merged$end, 30L)

  genes <- .make_gene_ranges(c(10L, 25L, 50L), c(30L, 40L, 60L))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 2L)
  expect_identical(merged$start, c(10L, 50L))
})

test_that(".merge_overlapping_gene_ranges handles edge cases", {
  # Empty input
  genes <- .make_gene_ranges(integer(0), integer(0))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 0L)

  # Non-overlapping genes
  genes <- .make_gene_ranges(c(10L, 50L, 100L), c(20L, 60L, 110L))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 3L)

  # Single gene
  genes <- .make_gene_ranges(100L, 200L)
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 1L)
})

test_that(".merge_overlapping_gene_ranges sorts results", {
  genes <- .make_gene_ranges(c(50L, 10L, 30L), c(60L, 20L, 40L))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(merged$start, sort(merged$start))
})

# .subdivide_oversized tests
test_that(".subdivide_oversized identifies needed splits", {
  # No splits needed
  result <- .subdivide_oversized(c(0L, 100L, 200L), limit = 250L)
  expect_identical(result, integer(0L))

  # Splits needed
  result <- .subdivide_oversized(c(0L, 500L), limit = 100L)
  expect_true(length(result) > 0L)
  expect_true(all(result > 0L & result < 500L))
})

test_that(".subdivide_oversized respects limits", {
  result <- .subdivide_oversized(c(0L, 500L, 1000L), limit = 150L)
  all_points <- sort(c(c(0L, 500L, 1000L), result))
  widths <- diff(all_points)
  expect_true(all(widths <= 150L))
})

test_that(".subdivide_oversized returns deduplicated sorted results", {
  result <- .subdivide_oversized(c(0L, 600L), limit = 100L)
  expect_identical(result, sort(unique(result)))
})

# .shift_single_boundary tests
test_that(".shift_single_boundary handles boundaries outside genes", {
  gene_ranges <- .make_gene_ranges(c(100L, 200L), c(150L, 250L))
  result <- .shift_single_boundary(50L, gene_ranges, 300L, 1L)
  expect_identical(result, 50L)
})

test_that(".shift_single_boundary shifts inside genes", {
  gene_ranges <- .make_gene_ranges(50L, 100L)
  result <- .shift_single_boundary(75L, gene_ranges, 300L, 1L)
  expect_gt(result, 100L)
  expect_lte(result, 300L)

  result <- .shift_single_boundary(75L, gene_ranges, 300L, 5L)
  expect_identical(result, 100L + 5L)
})

test_that(".shift_single_boundary respects bounds", {
  gene_ranges <- .make_gene_ranges(250L, 290L)
  result <- .shift_single_boundary(275L, gene_ranges, 300L, 1L)
  expect_lt(result, 300L)
})

# .shift_boundaries_past_genes tests
test_that(".shift_boundaries_past_genes handles edge cases", {
  gene_ranges <- .make_gene_ranges(50L, 100L)
  result <- .shift_boundaries_past_genes(integer(0L), gene_ranges, 300L, 1L)
  expect_identical(result, integer(0L))

  gene_ranges <- .make_gene_ranges(integer(0), integer(0))
  boundaries <- c(50L, 150L, 250L)
  result <- .shift_boundaries_past_genes(boundaries, gene_ranges, 300L, 1L)
  expect_identical(result, boundaries)
})

test_that(".shift_boundaries_past_genes deduplicates and filters", {
  # Test 1: deduplicates
  gene_ranges <- .make_gene_ranges(50L, 100L)
  boundaries <- c(75L, 75L, 200L)
  result <- .shift_boundaries_past_genes(boundaries, gene_ranges, 300L, 1L)
  expect_identical(length(result), length(unique(result)))

  # Test 2: filters boundaries (only interior boundaries kept)
  boundaries <- c(50L, 100L, 200L)
  gene_ranges <- .make_gene_ranges(integer(0), integer(0))
  result <- .shift_boundaries_past_genes(boundaries, gene_ranges, 300L, 1L)

  # All results should be valid interior boundaries
  expect_true(all(result > 0L & result < 300L))
})

# .assert_feasible_limit tests
test_that(".assert_feasible_limit accepts feasible configs", {
  genes <- data.frame(
    chr = character(0),
    start = integer(0),
    end = integer(0)
  )
  expect_invisible(.assert_feasible_limit(genes, limit = 50000L, shift_by = 1L))

  genes <- data.frame(
    chr = c("chr1", "chr1"),
    start = c(1000L, 5000L),
    end = c(2000L, 6000L),
    stringsAsFactors = FALSE
  )
  expect_invisible(.assert_feasible_limit(genes, limit = 2000L, shift_by = 1L))
})

test_that(".assert_feasible_limit rejects insufficient limits", {
  genes <- data.frame(
    chr = "chr1",
    start = 1000L,
    end = 6000L,
    stringsAsFactors = FALSE
  )
  expect_error(
    .assert_feasible_limit(genes, limit = 1000L, shift_by = 1L),
    regexp = "smaller than largest"
  )
})

test_that(".assert_feasible_limit accounts for shift_by", {
  genes <- data.frame(
    chr = "chr1",
    start = 1000L,
    end = 2000L,
    stringsAsFactors = FALSE
  )
  # Gene width is 1000, with shift_by=1 becomes 1001
  expect_error(
    .assert_feasible_limit(genes, limit = 1000L, shift_by = 1L),
    regexp = "smaller than"
  )
})

# .process_single_region tests
test_that(".process_single_region handles basic regions", {
  region <- data.frame(
    chr = "chr1",
    start = 0L,
    end = 1000L,
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    chr = character(0),
    start = integer(0),
    end = integer(0),
    stringsAsFactors = FALSE
  )
  result <- .process_single_region(region, genes, limit = 250L, shift_by = 1L)
  expect_s3_class(result, "data.frame")
  expect_true(all(result$end - result$start <= 250L))
})

test_that(".process_single_region handles centromeres", {
  region <- data.frame(
    chr = "chr1",
    start = 100L,
    end = 1000L,
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    chr = character(0),
    start = integer(0),
    end = integer(0),
    stringsAsFactors = FALSE
  )
  result <- .process_single_region(region, genes, limit = 250L, shift_by = 1L)
  expect_true(any(result$start == 0L))
  expect_true(any(result$start == 100L))
})

test_that(".process_single_region rejects invalid coords", {
  region <- data.frame(
    chr = "chr1",
    start = 100L,
    end = 100L,
    stringsAsFactors = FALSE
  )
  genes <- data.frame(
    chr = character(0),
    start = integer(0),
    end = integer(0),
    stringsAsFactors = FALSE
  )
  expect_error(
    .process_single_region(region, genes, limit = 250L, shift_by = 1L),
    regexp = "must be < end"
  )
})

# .read_bed_file tests
test_that(".read_bed_file reads valid files", {
  tmpfile <- tempfile(fileext = ".bed")
  writeLines(
    c("chr1\t0\t100", "chr1\t100\t200", "chr2\t0\t500"),
    tmpfile
  )
  on.exit(unlink(tmpfile))

  result <- .read_bed_file(tmpfile, min_cols = 3L)
  expect_s3_class(result, "data.frame")
  expect_identical(nrow(result), 3L)
  expect_identical(colnames(result)[1:3], c("chr", "start", "end"))
})

test_that(".read_bed_file validates structure", {
  tmpfile <- tempfile(fileext = ".bed")
  writeLines("chr1\t0", tmpfile)
  on.exit(unlink(tmpfile))

  expect_error(
    .read_bed_file(tmpfile, min_cols = 3L),
    regexp = "at least 3 columns"
  )
})

test_that(".read_bed_file requires file existence", {
  expect_error(
    .read_bed_file("/nonexistent/path.bed"),
    regexp = "not found"
  )
})

# .read_gtf_file tests
test_that(".read_gtf_file reads valid files", {
  tmpfile <- tempfile(fileext = ".gtf")
  writeLines(
    'chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
    tmpfile
  )
  on.exit(unlink(tmpfile))

  result <- .read_gtf_file(tmpfile)
  expect_identical(ncol(result), 9L)
  expect_identical(
    colnames(result),
    c(
      "chr",
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

test_that(".read_gtf_file skips comments and validates", {
  tmpfile <- tempfile(fileext = ".gtf")
  writeLines(
    c(
      "# This is a comment",
      'chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1";'
    ),
    tmpfile
  )
  on.exit(unlink(tmpfile))

  result <- .read_gtf_file(tmpfile)
  expect_identical(nrow(result), 1L)
})

test_that(".read_gtf_file validates structure", {
  tmpfile <- tempfile(fileext = ".gtf")
  writeLines("chr1\tensembl\tgene", tmpfile)
  on.exit(unlink(tmpfile))

  expect_error(
    .read_gtf_file(tmpfile),
    regexp = "at least 9 columns"
  )
})

# .write_bed_file tests
test_that(".write_bed_file writes correctly", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  df <- data.frame(
    chr = c("chr1", "chr2"),
    start = c(0L, 100L),
    end = c(50L, 150L),
    stringsAsFactors = FALSE
  )

  result <- .write_bed_file(df, tmpfile)
  expect_true(file.exists(tmpfile))

  content <- readLines(tmpfile)
  expect_length(content, 2L)
  expect_false(any(grepl("\"", content, fixed = TRUE)))
})
