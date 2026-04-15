.make_gtf_with_regions <- function(
  chrs,
  starts,
  ends,
  ids,
  region_chrs = NULL,
  region_starts = NULL,
  region_ends = NULL
) {
  gtf <- data.frame(
    Chr = chrs,
    start = starts,
    end = ends,
    unique_ID = ids,
    stringsAsFactors = FALSE
  )

  if (is.null(region_chrs)) {
    return(list(gtf = gtf, regions = NULL))
  }

  regions <- data.frame(
    Chr = region_chrs,
    RegionStart = region_starts,
    RegionEnd = region_ends,
    NEW = sprintf("%s:%s-%s", region_chrs, region_starts, region_ends),
    stringsAsFactors = FALSE
  )

  list(gtf = gtf, regions = regions)
}

# .format_region_name tests
test_that(".format_region_name formats correctly", {
  expect_identical(
    .format_region_name("chr1", 100000, 200000),
    "chr1:100000-200000"
  )
  expect_identical(.format_region_name("chrX", 0, 1000), "chrX:0-1000")
  expect_identical(.format_region_name("chr1", 1L, 2L), "chr1:1-2")
})

# .match_features_to_regions tests
test_that(".match_features_to_regions matches correctly", {
  data <- .make_gtf_with_regions(
    c("chr1", "chr1", "chr1"),
    c(100L, 200L, 500L),
    c(150L, 250L, 550L),
    c("ID1", "ID2", "ID3"),
    c("chr1", "chr1"),
    c(0L, 400L),
    c(300L, 600L)
  )

  result <- .match_features_to_regions(data$gtf, data$regions)
  expect_identical(nrow(result), 3L)
  expect_true(all(result$unique_ID %in% c("ID1", "ID2", "ID3")))
})

test_that(".match_features_to_regions excludes mismatched features", {
  # Feature crosses boundary
  data <- .make_gtf_with_regions(
    c("chr1", "chr1"),
    c(100L, 250L),
    c(150L, 350L),
    c("ID1", "ID2"),
    "chr1",
    0L,
    300L
  )

  result <- .match_features_to_regions(data$gtf, data$regions)
  expect_identical(nrow(result), 1L)
  expect_identical(result$unique_ID, "ID1")

  # Feature before region
  data <- .make_gtf_with_regions(
    "chr1",
    100L,
    150L,
    "ID1",
    "chr1",
    200L,
    500L
  )

  result <- .match_features_to_regions(data$gtf, data$regions)
  expect_identical(nrow(result), 0L)
})

# .check_no_ambiguous_matches tests
test_that(".check_no_ambiguous_matches validates uniqueness", {
  matched <- data.frame(
    unique_ID = c("ID1", "ID2", "ID3"),
    stringsAsFactors = FALSE
  )
  expect_invisible(.check_no_ambiguous_matches(matched))

  matched <- data.frame(
    unique_ID = character(0)
  )
  expect_invisible(.check_no_ambiguous_matches(matched))
})

test_that(".check_no_ambiguous_matches detects duplicates", {
  matched <- data.frame(
    unique_ID = c("ID1", "ID1", "ID2"),
    stringsAsFactors = FALSE
  )
  expect_error(
    .check_no_ambiguous_matches(matched),
    regexp = "Ambiguous"
  )

  matched <- data.frame(
    unique_ID = c("ID1", "ID1", "ID2", "ID2"),
    stringsAsFactors = FALSE
  )
  expect_error(
    .check_no_ambiguous_matches(matched),
    regexp = "ID1.*ID2"
  )
})

# .check_all_features_matched tests
test_that(".check_all_features_matched accepts complete matches", {
  matched <- data.frame(unique_ID = c("ID1", "ID2"), stringsAsFactors = FALSE)
  gtf <- data.frame(unique_ID = c("ID1", "ID2"), stringsAsFactors = FALSE)

  expect_invisible(.check_all_features_matched(matched, gtf))
})

test_that(".check_all_features_matched detects missing features", {
  matched <- data.frame(unique_ID = "ID1", stringsAsFactors = FALSE)
  gtf <- data.frame(
    unique_ID = c("ID1", "ID2", "ID3"),
    Chr = c("chr1", "chr1", "chr1"),
    feature = c("gene", "gene", "gene"),
    start = c(0L, 1000L, 2000L),
    end = c(100L, 2000L, 3000L),
    stringsAsFactors = FALSE
  )

  expect_error(
    .check_all_features_matched(matched, gtf),
    regexp = "do not fit within"
  )

  # Checks count and examples
  expect_error(
    .check_all_features_matched(matched, gtf),
    regexp = "2.*Examples"
  )
})

# .build_output_gtf tests
test_that(".build_output_gtf reorders and selects columns", {
  matched <- data.frame(
    unique_ID = c("ID2", "ID1"),
    NEW = c("chr1:0-300", "chr1:0-200"),
    source = c("source", "source"),
    feature = c("gene", "gene"),
    start = c(200L, 100L),
    end = c(300L, 150L),
    score = c(".", "."),
    strand = c("+", "+"),
    frame = c(".", "."),
    attribute = c('gene_id "G2";', 'gene_id "G1";'),
    stringsAsFactors = FALSE
  )

  gtf_file <- data.frame(
    unique_ID = c("ID1", "ID2"),
    stringsAsFactors = FALSE
  )

  result <- .build_output_gtf(matched, gtf_file)
  expect_identical(result$NEW, c("chr1:0-200", "chr1:0-300"))
  expect_identical(colnames(result)[1], "Chr")
  expect_identical(result$Chr, c("chr1:0-200", "chr1:0-300"))
})

test_that(".build_output_gtf reorders and selects columns", {
  matched <- data.frame(
    unique_ID = c("ID2", "ID1"),
    NEW = c("chr1:0-300", "chr1:0-200"),
    source = c("source", "source"),
    feature = c("gene", "gene"),
    start = c(200L, 100L),
    end = c(300L, 150L),
    score = c(".", "."),
    strand = c("+", "+"),
    frame = c(".", "."),
    attribute = c('gene_id "G2";', 'gene_id "G1";'),
    stringsAsFactors = FALSE
  )

  gtf_file <- data.frame(
    unique_ID = c("ID1", "ID2"),
    stringsAsFactors = FALSE
  )

  result <- .build_output_gtf(matched, gtf_file)
  # Chr column should have sorted NEW values
  expect_identical(result$Chr, c("chr1:0-200", "chr1:0-300"))
  expect_identical(colnames(result)[1L], "Chr")
})
