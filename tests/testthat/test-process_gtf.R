# Test helper function
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
    source = "source",
    feature = "gene",
    score = ".",
    strand = "+",
    frame = ".",
    attribute = 'gene_id "test";',
    stringsAsFactors = FALSE
  )

  if (is.null(region_chrs)) {
    return(list(gtf = gtf, regions = NULL))
  }

  regions <- data.frame(
    Chr = region_chrs,
    RegionStart = region_starts,
    RegionEnd = region_ends,
    stringsAsFactors = FALSE
  )

  list(gtf = gtf, regions = regions)
}

# .assign_coords tests
test_that(".assign_coords adds NEW column correctly", {
  df <- data.frame(
    Chr = c("1H", "2H"),
    start = c(100, 200),
    end = c(150, 250),
    stringsAsFactors = FALSE
  )

  result <- .assign_coords(df, c("1H", "2H"), c(0, 100), c(200, 300))
  expect_true("NEW" %in% colnames(result))
  expect_identical(result$NEW, c("1H-0-200", "2H-100-300"))
})

test_that(".assign_coords handles empty dataframes", {
  df <- data.frame(
    Chr = character(0),
    start = numeric(0),
    end = numeric(0),
    stringsAsFactors = FALSE
  )

  result <- .assign_coords(df, character(0), numeric(0), numeric(0))
  expect_identical(nrow(result), 0L)
  expect_true("NEW" %in% colnames(result) || nrow(result) == 0L)
})

test_that(".assign_coords preserves input rows", {
  df <- data.frame(
    Chr = c("chr1", "chr1", "chr1"),
    start = c(100, 200, 300),
    end = c(150, 250, 350),
    stringsAsFactors = FALSE
  )

  result <- .assign_coords(df, df$Chr, c(0, 50, 100), c(200, 250, 350))
  expect_identical(nrow(result), 3L)
  expect_identical(result$start, df$start)
  expect_identical(result$end, df$end)
})

# .recalculate_above_centromere_coords tests
test_that(".recalculate_above_centromere_coords recalculates coordinates", {
  df <- data.frame(
    Chr = "1H",
    start = 100L,
    end = 150L,
    CentrStart = 50L,
    ChrEND = 200L,
    stringsAsFactors = FALSE
  )

  result <- .recalculate_above_centromere_coords(df)

  # Expected calculation:
  # feature_length = 150 - 100 = 50
  # dist_chr_end_feat_end = 200 - 150 = 50
  # dist_centr_chr_end = 200 - 50 = 150
  # new_start = 150 - 50 - 50 = 50
  # new_end = 150 - 50 = 100
  expect_identical(result$start, 50L)
  expect_identical(result$end, 100L)
})

test_that(".recalculate_above_centromere_coords handles multiple rows", {
  df <- data.frame(
    Chr = c("1H", "1H"),
    start = c(100L, 110L),
    end = c(150L, 160L),
    CentrStart = c(50L, 50L),
    ChrEND = c(200L, 200L),
    stringsAsFactors = FALSE
  )

  result <- .recalculate_above_centromere_coords(df)
  expect_identical(nrow(result), 2L)
  expect_true(all(result$start >= 0L))
  expect_true(all(result$end >= result$start))
})

test_that(".recalculate_above_centromere_coords handles edge case", {
  # Feature at the very end
  df <- data.frame(
    Chr = "1H",
    start = 150L,
    end = 200L,
    CentrStart = 50L,
    ChrEND = 200L,
    stringsAsFactors = FALSE
  )

  result <- .recalculate_above_centromere_coords(df)
  # feature_length = 200 - 150 = 50
  # dist_chr_end_feat_end = 200 - 200 = 0
  # dist_centr_chr_end = 200 - 50 = 150
  # new_start = 150 - 0 - 50 = 100
  # new_end = 150 - 0 = 150
  expect_identical(result$start, 100L)
  expect_identical(result$end, 150L)
})

# process_gtf integration tests
test_that("process_gtf processes basic GTF correctly", {
  skip_if_not_installed("cli")

  # Create temporary files
  gtf_content <- "1H\tsource\tgene\t0\t20\t.\t+\t.\tgene_id \"gene1\";
1H\tsource\tgene\t100\t120\t.\t+\t.\tgene_id \"gene2\";"

  bed_content <- "1H\t50\t150"

  gtf_file <- tempfile(fileext = ".gtf")
  bed_file <- tempfile(fileext = ".bed")

  writeLines(gtf_content, gtf_file)
  writeLines(bed_content, bed_file)

  on.exit({
    unlink(gtf_file)
    unlink(bed_file)
  })

  result <- process_gtf(
    gtf_path = gtf_file,
    centromere_path = bed_file,
    unchar_region = "CAJHDD.*",
    mito = "Mt",
    chloro = "Pt"
  )

  expect_true(is.data.frame(result))
  expect_true("NEW" %in% colnames(result))
  expect_identical(nrow(result), 2L)
})

test_that("process_gtf handles special contigs", {
  skip_if_not_installed("cli")

  gtf_content <- "Mt\tsource\tgene\t1\t10\t.\t+\t.\tgene_id \"MT1\";
Pt\tsource\tgene\t2\t8\t.\t+\t.\tgene_id \"PT1\";"

  bed_content <- "1H\t50\t150"

  gtf_file <- tempfile(fileext = ".gtf")
  bed_file <- tempfile(fileext = ".bed")

  writeLines(gtf_content, gtf_file)
  writeLines(bed_content, bed_file)

  on.exit({
    unlink(gtf_file)
    unlink(bed_file)
  })

  result <- process_gtf(
    gtf_path = gtf_file,
    centromere_path = bed_file,
    unchar_region = "CAJHDD.*",
    mito = "Mt",
    chloro = "Pt"
  )

  expect_true(is.data.frame(result))
  expect_identical(nrow(result), 2L)
  expect_true(any(grepl("Mt-", result$NEW)))
  expect_true(any(grepl("Pt-", result$NEW)))
})

test_that("process_gtf returns correct columns", {
  skip_if_not_installed("cli")

  gtf_content <- "1H\tsource\tgene\t0\t20\t.\t+\t.\tgene_id \"gene1\";"
  bed_content <- "1H\t50\t150"

  gtf_file <- tempfile(fileext = ".gtf")
  bed_file <- tempfile(fileext = ".bed")

  writeLines(gtf_content, gtf_file)
  writeLines(bed_content, bed_file)

  on.exit({
    unlink(gtf_file)
    unlink(bed_file)
  })

  result <- process_gtf(
    gtf_path = gtf_file,
    centromere_path = bed_file,
    unchar_region = "CAJHDD.*",
    mito = "Mt",
    chloro = "Pt"
  )

  expected_cols <- c(
    "NEW",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )
  expect_identical(colnames(result), expected_cols)
})
