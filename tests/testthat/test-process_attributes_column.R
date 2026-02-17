test_that("process_attributes_column() edits 9th column in GTF file", {
  
  # output dir:
  out_path2 <- test_path("testdata/testoutput/") 
  dir.create(out_path2, showWarnings = FALSE)
  
  # input param for process_attributes_column()
  keep_attrib= c('gene_id', 'transcript_id', 'gene_name')
  genome_n ="BarleyMOREXV3_toy"
  genome_ver ="v51_ENSEMBL"

 #  Step 1: Run process_gtf() function 
 result1 <- process_gtf(
    gtf_path = test_path("testdata/A3_toy_all_scenarios_2chr.gtf"),
    centromere_path = test_path("testdata/testoutput/RECAL1_FINAL_SPLIT_coord_for_GTF.bed"),
    unchar_region = "^CAJHDD",   #### "^ChrUn"
    mito = "^Mt$",               #### "^ChrM$" 
    chloro = "^Pt$"              #### "^ChrC$" "^NA$" # if geenome has no chloroplast(mammalian) )
    )

  #  Step 2: Run process_attributes_column() function 
  expect_no_error(
    result2 <- process_attributes_column(
    result = result1, 
    keep_attributes = keep_attrib, 
    out_path = out_path2, 
    genome_name = genome_n, 
    genome_version = genome_ver
      )
    )
  # Optional expectations 
  expect_true(!is.null(result2))
  })

