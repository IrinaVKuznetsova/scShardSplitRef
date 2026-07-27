# Tests for messaging.R (.scs_get_verbosity / .scs_emit_info),
# scShardSplitRef-options.R (scShardSplitRef_options()), and the verbosity
# mapping in zzz.R (.map_verbosity()).
#
# None of these were previously covered by any test file, even though
# verbosity gates every cli_inform()/cli_alert_info() call in
# determine_split_regions() and process_gtf().

# .scs_get_verbosity ----------------------------------------------------------

test_that(".scs_get_verbosity returns the option value when it is valid", {
  withr::local_options(scShardSplitRef.verbosity = "quiet")
  expect_identical(.scs_get_verbosity(), "quiet")

  withr::local_options(scShardSplitRef.verbosity = "minimal")
  expect_identical(.scs_get_verbosity(), "minimal")

  withr::local_options(scShardSplitRef.verbosity = "verbose")
  expect_identical(.scs_get_verbosity(), "verbose")
})

test_that(".scs_get_verbosity defaults to 'verbose' when the option is unset", {
  withr::local_options(scShardSplitRef.verbosity = NULL)
  expect_identical(.scs_get_verbosity(), "verbose")
})

test_that(".scs_get_verbosity falls back to 'verbose' for invalid values", {
  withr::local_options(scShardSplitRef.verbosity = "loud")
  expect_identical(.scs_get_verbosity(), "verbose")

  withr::local_options(scShardSplitRef.verbosity = 1L)
  expect_identical(.scs_get_verbosity(), "verbose")
})

# .scs_emit_info ----------------------------------------------------------------

test_that(".scs_emit_info() is TRUE only when verbosity is 'verbose'", {
  withr::local_options(scShardSplitRef.verbosity = "verbose")
  expect_true(.scs_emit_info())

  withr::local_options(scShardSplitRef.verbosity = "minimal")
  expect_false(.scs_emit_info())

  withr::local_options(scShardSplitRef.verbosity = "quiet")
  expect_false(.scs_emit_info())
})

# scShardSplitRef_options() -----------------------------------------------------

test_that("scShardSplitRef_options() with no args returns only scShardSplitRef.* options", {
  withr::local_options(scShardSplitRef.verbosity = "minimal")

  result <- scShardSplitRef_options()

  expect_type(result, "list")
  expect_true("scShardSplitRef.verbosity" %in% names(result))
  expect_identical(result$scShardSplitRef.verbosity, "minimal")
  expect_true(all(startsWith(names(result), "scShardSplitRef.")))
})

test_that("scShardSplitRef_options(...) sets options that getOption() then reflects", {
  old <- getOption("scShardSplitRef.verbosity")
  on.exit(options(scShardSplitRef.verbosity = old), add = TRUE)

  scShardSplitRef_options(scShardSplitRef.verbosity = "quiet")
  expect_identical(getOption("scShardSplitRef.verbosity"), "quiet")

  scShardSplitRef_options(scShardSplitRef.verbosity = "verbose")
  expect_identical(getOption("scShardSplitRef.verbosity"), "verbose")
})

# .map_verbosity ------------------------------------------------------------

test_that(".map_verbosity maps 'quiet' to fully silent derived options", {
  m <- .map_verbosity("quiet")
  expect_identical(m$rlib_message_verbosity, "quiet")
  expect_identical(m$rlib_warning_verbosity, "quiet")
  expect_identical(m$warn, -1L)
  expect_false(m$datatable.showProgress)
})

test_that(".map_verbosity maps 'minimal' to messages-off, warnings-on", {
  m <- .map_verbosity("minimal")
  expect_identical(m$rlib_message_verbosity, "minimal")
  expect_identical(m$rlib_warning_verbosity, "verbose")
  expect_identical(m$warn, 0L)
  expect_false(m$datatable.showProgress)
})

test_that(".map_verbosity maps 'verbose' to everything on", {
  m <- .map_verbosity("verbose")
  expect_identical(m$rlib_message_verbosity, "verbose")
  expect_identical(m$rlib_warning_verbosity, "verbose")
  expect_identical(m$warn, 0L)
  expect_true(m$datatable.showProgress)
})

test_that(".map_verbosity treats NULL and unrecognized values as 'verbose'", {
  expect_identical(.map_verbosity(NULL), .map_verbosity("verbose"))
  expect_identical(.map_verbosity("nonsense"), .map_verbosity("verbose"))
})

# End-to-end: verbosity actually silences/emits pipeline messages --------------

test_that("determine_split_regions() emits messages only when verbose", {
  bed <- withr::local_tempfile(fileext = ".bed")
  gtf <- withr::local_tempfile(fileext = ".gtf")
  out <- withr::local_tempfile(fileext = ".bed")

  writeLines("chr1\t0\t1000", bed)
  writeLines('chr1\tsrc\tgene\t500\t600\t.\t+\t.\tgene_id "G1";', gtf)

  withr::local_options(scShardSplitRef.verbosity = "verbose")
  expect_message(
    determine_split_regions(
      bed = bed,
      gtf = gtf,
      output_bed = out,
      limit = 2000L
    ),
    "Preparing split regions"
  )

  withr::local_options(scShardSplitRef.verbosity = "minimal")
  expect_no_message(
    determine_split_regions(
      bed = bed,
      gtf = gtf,
      output_bed = out,
      limit = 2000L
    )
  )

  withr::local_options(scShardSplitRef.verbosity = "quiet")
  expect_no_message(
    determine_split_regions(
      bed = bed,
      gtf = gtf,
      output_bed = out,
      limit = 2000L
    )
  )
})

test_that("process_gtf() emits messages only when verbose", {
  reg_file <- withr::local_tempfile(fileext = ".bed")
  gtf_file <- withr::local_tempfile(fileext = ".gtf")
  out_path <- withr::local_tempdir()

  writeLines("chr1\t0\t100", reg_file)
  writeLines('chr1\tsrc\tgene\t10\t20\t.\t+\t.\tgene_id "G1";', gtf_file)

  withr::local_options(scShardSplitRef.verbosity = "verbose")
  expect_message(
    process_gtf(
      split_regions_bed = reg_file,
      gtf = gtf_file,
      genome_name = "verbosetest",
      genome_version = "v1",
      out_path = out_path
    ),
    "Validating inputs"
  )

  withr::local_options(scShardSplitRef.verbosity = "quiet")
  expect_no_message(
    process_gtf(
      split_regions_bed = reg_file,
      gtf = gtf_file,
      genome_name = "quiettest",
      genome_version = "v1",
      out_path = out_path
    )
  )
})
