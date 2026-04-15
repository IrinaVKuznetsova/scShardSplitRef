# process_attributes_column integration tests
test_that("process_attributes_column filters correctly", {
  tmpdir <- tempdir()

  result_df <- data.frame(
    attribute = c(
      'gene_id "G1"; transcript_id "T1"; gene_name "FOO"',
      'gene_id "G2"; transcript_id "T2"; gene_name "BAR"'
    ),
    stringsAsFactors = FALSE
  )

  process_attributes_column(
    result = result_df,
    keep_attributes = c("gene_id", "transcript_id"),
    out_path = tmpdir,
    genome_name = "test",
    genome_version = "v1"
  )

  output_file <- file.path(tmpdir, "B1_FINAL_MODIFIED_GTF_test_v1.gtf")
  expect_true(file.exists(output_file))
  on.exit(unlink(output_file))

  content <- readLines(output_file)
  expect_true(all(grepl("gene_id", content)))
  expect_false(any(grepl("gene_name", content)))
})

test_that("process_attributes_column validates inputs", {
  result_df <- data.frame(
    attribute = 'gene_id "G1";',
    stringsAsFactors = FALSE
  )

  expect_error(
    process_attributes_column(
      result = result_df,
      keep_attributes = "gene_id",
      out_path = "/nonexistent/path",
      genome_name = "test",
      genome_version = "v1"
    ),
    regexp = "does not exist"
  )
})

test_that("process_attributes_column warns on incomplete filtering", {
  tmpdir <- tempdir()

  result_df <- data.frame(
    attribute = c('gene_id "G1";', 'transcript_id "T1";'),
    stringsAsFactors = FALSE
  )

  expect_warning(
    process_attributes_column(
      result = result_df,
      keep_attributes = "gene_id",
      out_path = tmpdir,
      genome_name = "test",
      genome_version = "v1"
    ),
    regexp = "no matching attributes"
  )

  output_file <- file.path(tmpdir, "B1_FINAL_MODIFIED_GTF_test_v1.gtf")
  on.exit(unlink(output_file))
})

test_that(".filter_attributes keeps specified attributes", {
  attributes <- c(
    'gene_id "G1"; gene_name "GENE1"; gene_type "protein_coding"',
    'gene_id "G2"; gene_name "GENE2"; gene_type "lncRNA"'
  )

  result <- .filter_attributes(attributes, c("gene_id", "gene_name"))

  expect_true(all(grepl("gene_id", result)))
  expect_true(all(grepl("gene_name", result)))
  expect_false(any(grepl("gene_type", result)))
})

test_that(".check_filtered_completeness warns on empty results", {
  filtered <- c("", "gene_id \"G2\"")
  original <- c("gene_id \"G1\"", "gene_id \"G2\"")

  expect_warning(
    .check_filtered_completeness(filtered, original),
    regexp = "matching"
  )
})

test_that(".build_output_filename creates correct path", {
  result <- .build_output_filename(".", "hg38", "v1")
  expect_true(grepl("B1_FINAL_MODIFIED_GTF_hg38_v1.gtf", result))
})

test_that(".write_filtered_gtf writes file", {
  tmpfile <- tempfile(fileext = ".gtf")
  on.exit(unlink(tmpfile))

  df <- data.frame(
    chr = "chr1",
    start = 1000L,
    end = 2000L,
    stringsAsFactors = FALSE
  )

  .write_filtered_gtf(df, tmpfile)
  expect_true(file.exists(tmpfile))
})
