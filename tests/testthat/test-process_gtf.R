test_that("process_gtf() modifies GTF file accounting for centromere/split coordinates", {
 
  expect_no_error({

    result <- process_gtf(
    gtf_path = test_path("testdata/A3_toy_all_scenarios_2chr.gtf"),
    centromere_path = test_path("testdata/testoutput/RECAL1_FINAL_SPLIT_coord_for_GTF.bed"),
    unchar_region = "^CAJHDD",   #### "^ChrUn"
    mito = "^Mt$",               #### "^ChrM$" 
    chloro = "^Pt$"              #### "^ChrC$" "^NA$" # if geenome has no chloroplast(mammalian) )
      )
    })
    expect_true(!is.null(result))  
})