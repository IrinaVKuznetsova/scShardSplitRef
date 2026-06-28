test_that(".build_boundaries respects size limit", {
  expect_identical(.build_boundaries(0L, 100L, 200L), integer(0L))
  result <- .build_boundaries(0L, 1000L, 250L)
  expect_type(result, "integer")
  expect_gt(length(result), 0L)
})

test_that(".build_boundaries creates valid chunks", {
  boundaries <- .build_boundaries(0L, 1000L, 250L)
  points <- c(0L, boundaries, 1000L)
  widths <- diff(points)
  expect_true(all(widths <= 250L))
})

test_that(".build_boundaries handles exact multiples", {
  boundaries <- .build_boundaries(0L, 400L, 100L)
  points <- c(0L, boundaries, 400L)
  expect_identical(sum(diff(points)), 400L)
})

test_that(".validate_bed_coordinates accepts valid coords", {
  expect_invisible(.validate_bed_coordinates(0L, 100L, "chr1"))
  expect_invisible(.validate_bed_coordinates(1000L, 2000L, "chr1"))
})

test_that(".validate_bed_coordinates rejects invalid coords", {
  expect_error(
    .validate_bed_coordinates(100L, 100L, "chr1"),
    regexp = "Invalid BED region"
  )
})

test_that(".boundaries_to_df creates correct structure", {
  df <- .boundaries_to_df("chr1", 0L, 100L, c(30L, 60L))
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), 3L)
  expect_identical(colnames(df), c("chr", "start", "end"))
  expect_identical(df$start, c(0L, 30L, 60L))
  expect_identical(df$end, c(30L, 60L, 100L))
})

test_that(".boundaries_to_df handles edge cases", {
  df <- .boundaries_to_df("chr1", 0L, 50L, integer(0L))
  expect_identical(nrow(df), 1L)

  df <- .boundaries_to_df("chrX", 0L, 50L, 25L)
  expect_identical(df$chr, c("chrX", "chrX"))
})

test_that(".boundaries_to_df ensures non-overlapping intervals", {
  df <- .boundaries_to_df("chr1", 0L, 1000L, c(250L, 500L, 750L))
  for (i in seq_len(nrow(df) - 1)) {
    expect_identical(df$end[i], df$start[i + 1])
  }
})

test_that(".split_region_core respects size limit", {
  result <- .split_region_core("chr1", 0L, 1000L, 250L)
  expect_s3_class(result, "data.frame")
  expect_true(all(result$end - result$start <= 250L))

  result <- .split_region_core("chr1", 0L, 100L, 500L)
  expect_identical(nrow(result), 1L)
})

test_that(".split_region_core preserves chromosome", {
  result <- .split_region_core("chrY", 0L, 100L, 50L)
  expect_identical(unique(result$chr), "chrY")
})

test_that(".merge_overlapping_gene_ranges merges overlaps", {
  genes <- data.frame(start = c(10L, 15L), end = c(30L, 25L))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 1L)
  expect_identical(merged$start, 10L)
})

test_that(".merge_overlapping_gene_ranges handles edge cases", {
  genes <- data.frame(start = integer(0), end = integer(0))
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 0L)
})

test_that(".subdivide_oversized identifies needed splits", {
  # No splits needed - all intervals within limit
  result <- .subdivide_oversized(c(0L, 100L, 200L), limit = 250L)
  expect_identical(result, integer(0))

  # Splits needed - 500bp interval with 100bp limit
  # Should create ~4 boundaries to split into 5 chunks of ~100bp each
  result <- .subdivide_oversized(c(0L, 500L), limit = 100L)

  expect_gt(length(result), 0L)
  expect_true(all(result > 0L & result < 500L))

  # Verify resulting intervals respect limit
  points <- c(0L, sort(result), 500L)
  widths <- diff(points)
  expect_true(all(widths <= 100L))
})

test_that(".shift_single_boundary handles boundaries outside genes", {
  gene_ranges <- data.frame(start = c(100L, 200L), end = c(150L, 250L))
  result <- .shift_single_boundary(50L, gene_ranges, 300L, 1L, 1L)
  expect_identical(result, 50L)
})

test_that(".shift_boundaries_past_genes handles edge cases", {
  gene_ranges <- data.frame(start = 50L, end = 100L)
  result <- .shift_boundaries_past_genes(integer(0L), gene_ranges, 300L, 1L)
  expect_identical(result, integer(0L))
})

test_that(".assert_feasible_limit accepts feasible configs", {
  genes <- data.frame(
    chr = character(0),
    start = integer(0),
    end = integer(0),
    stringsAsFactors = FALSE
  )
  expect_invisible(.assert_feasible_limit(
    genes,
    limit = 50000L,
    clearance = 1L
  ))
})

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
})

test_that(".read_bed_file reads valid files", {
  tmpfile <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t0\t100", "chr1\t100\t200"), tmpfile)
  on.exit(unlink(tmpfile))

  result <- .read_bed_file(tmpfile, min_cols = 3L)
  expect_s3_class(result, "data.frame")
  expect_identical(nrow(result), 2L)
})

test_that(".read_gtf_file reads valid files", {
  tmpfile <- tempfile(fileext = ".gtf")
  writeLines('chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "G1";', tmpfile)
  on.exit(unlink(tmpfile))

  result <- .read_gtf_file(tmpfile)
  expect_identical(ncol(result), 9L)
})

test_that(".write_bed_file writes BED file and returns invisibly", {
  tmpfile <- tempfile(fileext = ".bed")
  on.exit(unlink(tmpfile))

  df <- data.frame(
    chr = c("chr1", "chr2"),
    start = c(0L, 100L),
    end = c(50L, 150L)
  )

  result <- .write_bed_file(df, tmpfile)

  expect_true(file.exists(tmpfile))
  expect_true(result)

  # Verify file contents
  content <- readLines(tmpfile)
  expect_length(content, 2L)
  expect_false(any(grepl("\"", content, fixed = TRUE)))
})
