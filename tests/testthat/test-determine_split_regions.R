test_that(".build_boundaries respects limit", {
  expect_identical(.build_boundaries(0L, 100L, 200L), integer(0L))
  expect_gt(.build_boundaries(0L, 100L, 30L), 0L)
})

test_that(".split_region_core rejects invalid BED", {
  bad <- data.frame(
    chr = "chr1",
    start = 10,
    end = 10,
    stringsAsFactors = FALSE
  )
  expect_error(.split_region_core(bad, 100))
})

test_that(".split_region_core produces valid partitions", {
  r <- data.frame(
    chr = "chr1",
    start = 0L,
    end = 100L,
    stringsAsFactors = FALSE
  )
  out <- .split_region_core(r, limit = 30L)
  expect_true(all(out$end > out$start))
  expect_lte(out$end - out$start, 30L)
})

test_that(".merge_overlapping_gene_ranges merges overlaps", {
  genes <- data.frame(
    chr = "chr1",
    start = c(10, 15, 30),
    end = c(20, 25, 40)
  )
  merged <- .merge_overlapping_gene_ranges(genes)
  expect_identical(nrow(merged), 2L)
})
