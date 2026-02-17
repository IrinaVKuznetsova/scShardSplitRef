test_that("prepare_centromere_split() handles gene-centromere overlaps", {

  out_path <- test_path("testdata/testoutput/") 
  dir.create(out_path, showWarnings = FALSE)
  
  expect_no_error(
    prepare_centromere_split(
      gtf_path = "testdata/A3_toy_all_scenarios_2chr.gtf",
      centromere_path = "testdata/IN0_toy_centromeres_for_gtf.bed",
      out_path = out_path,
      contig_prefix = "CAJHDD", 
      mt_name = "Mt", 
      pt_name = "Pt"
      )
    )
  })