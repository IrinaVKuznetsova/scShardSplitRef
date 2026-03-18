#' Check for gene feature overlaps with centromere coordinates and prepare new split coordinates if needed
#'
#' Converts a BED file to a format suitable for splitting GTF files
#'
#' @param df A data frame containing centromere coordinates with columns "Chr", "CentrStart", and "ChrEND".
#' @param contig_prefix A string prefix to identify contigs in the chromosome names (
#' @param mt_name A string representing the name of the mitochondrial chromosome (e.g., "Mt").
#' @param pt_name A string representing the name of the plastid chromosome (e.g., "Pt").
#'
#' @import dplyr
#' @import stringr
#' @importFrom utils read.table write.table
#'
#' @examples
#' prepare_centromere_BED(df,
#'                        contig_prefix = "CAJHDD",
#'                        mt_name = "Mt",
#'                        pt_name = "Pt")
#'
#' @author Irina Kuznetsova, \email{irina.Kuznetsova@@curtin.edu.au}
#' @export
#'
#################################################################################################################################################
# FUNCTION-1 checks for the rare case where a gene feature lies on the centromere coordinate and recalculates the split coordinate.
# It takes the end coordinate of any gene that overlaps the centromere and shifts it 10 bp downstream to establish a new split coordinate.
#
# FUNCTION-0 prepares centromere coordinates in the correct format for FASTA splitting.
#################################################################################################################################################

prepare_centromere_BED = function(df, contig_prefix, mt_name, pt_name) {
  df_result = df %>%
    rowwise() %>%
    do({
      chr = .$Chr
      start = .$CentrStart
      end = .$ChrEND

      # Identify contigs and Mt and Pt
      contig = str_detect(chr, paste0("^", contig_prefix))
      organelle = chr %in% c(mt_name, pt_name)
      # check if either condition is TRUE
      if (contig == TRUE || organelle == TRUE) {
        # Keep unchanged
        data.frame(Chr = chr, CentrStart = start, ChrEND = end)
      } else {
        # Split only if CentrStart > 0
        if (start > 0) {
          data.frame(
            Chr = chr,
            CentrStart = c(0, start),
            ChrEND = c(start, end)
          )
        } else {
          data.frame(Chr = chr, CentrStart = start, ChrEND = end)
        }
      }
    }) %>%
    ungroup()
  return(df_result)
}


#' Check if a gene feature overlaps with centromere coordinates and prepare new split coordinates if needed
#'
#' Check if GTF file has overlapping gene features with centromere coordinates
#' and prepare new split coordinates if that's the case.
#'
#' @param gtf_path A string representing the file path to the GTF file.
#' @param centromere_path A string representing the file path to the centromere
#'   coordinates BED file.
#' @param out_path A string representing the directory path where the output
#'  files will be saved.
#' @param mt_name A string representing the name of the mitochondrial chromosome
#'  (*e.g.*, "Mt").
#' @param pt_name A string representing the name of the plastid chromosome
#'  (*e.g.*, "Pt").
#'
#' @import dplyr
#' @import stringr
#' @importFrom utils read.table write.table
#'
#' @examples
#' prep_centromere_split(gtf_path = "/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/A0_toy_all_scenarios_2chr.gtf,   #A0_toy_all_scenarios.gtf",     ###/A0_toy_all_scenarios.gtf",
#'                       centromere_path = "/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/A0_toy_centromeres_for_gtf.bed",
#'                       out_path = "/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/",
#'                       contig_prefix = "CAJHDD",
#'                       mt_name = "Mt",
#'                       pt_name = "Pt")
#'
#' @author Irina Kuznetsova, \email{irina.Kuznetsova@@curtin.edu.au}
#' @export
#'
prepare_centromere_split = function(
  gtf_path,
  centromere_path,
  out_path,
  contig_prefix,
  mt_name,
  pt_name
) {
  # --------------------------------------------------
  # 1 Read GTF file, add column names
  # --------------------------------------------------
  gtf_file = read.table(file = gtf_path, header = FALSE, sep = "\t")
  colnames(gtf_file) = c(
    "Chr",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )
  gtf_file$unique_ID = paste0("ID", seq_len(nrow(gtf_file)))

  # --------------------------------------------------
  # 2 Read centromere positions file (BED format)
  # --------------------------------------------------
  read_cent = read.table(file = centromere_path, sep = "\t")
  colnames(read_cent) = c("Chr", "CentrStart", "ChrEND")

  # Convert original BED to suit FASTA splitting
  orig_For_FASTA = prep_centromere_BED(
    read_cent,
    contig_prefix,
    mt_name,
    pt_name
  )
  write.table(
    orig_For_FASTA,
    file = paste0(out_path, "IN0_original_centromeres_for_FASTA.bed"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  # ---------------------------------------------------------------------------
  # 3 Combine two data frames
  # --------------------------------------------------------------------------
  df_merged_chr <- merge(gtf_file, read_cent, by = "Chr")

  # --------------------------------------------------------------------------
  # 4 Check for overlapping gene features
  # --------------------------------------------------------------------------
  check_genefeat_centr_overlap = df_merged_chr %>%
    dplyr::filter(start < CentrStart & end > CentrStart)

  # if (nrow(check_genefeat_centr_overlap) > 0){
  #   print(c(paste("g_start:", check_genefeat_centr_overlap$start,
  #           "g_end:", check_genefeat_centr_overlap$end,
  #           "chr_Start:", check_genefeat_centr_overlap$CentrStart,
  #           "chr_End:", check_genefeat_centr_overlap$ChrEND, sep = " "),
  #           "Attention this genome has gene features coordinates overlaping with the centromere cooridnate in",
  #           check_genefeat_centr_overlap$Chr[1]    ))
  # }

  # Prepare empty object to store new split coordinates
  BEDnew_split_coord_for_overlap_genes_with_centrom <- data.frame()

  # --------------------------------------------------------------------------
  # 5 Handle overlaps
  # --------------------------------------------------------------------------
  if (nrow(check_genefeat_centr_overlap) > 0) {
    message(
      "\nOverlapping gene features with centromere coordinates detected!\n",
      "re-calculating split coordinates for these chromosomes\n",
      "Use:\n",
      out_path,
      "RECAL1_FINAL_SPLIT_coord_for_GTF.bed as the input to modify GTF file.\n",
      "and:\n",
      out_path,
      "RECAL1_FINAL_SPLIT_coord_for_FASTA.bed",
      "\n as the input for splitting FASTA."
    )

    # 5.1 Shift centromere split coordinate 10 bp downstream of the gene feature end coordinate
    new_split_coord_for_overlap_genes_with_centrom = check_genefeat_centr_overlap$end +
      10

    # 5.2 Create new data frame with updated split coordinates
    new_split_coord_for_overlap_genes_with_centrom = data.frame(
      Chr = check_genefeat_centr_overlap$Chr,
      source = check_genefeat_centr_overlap$source,
      feature = check_genefeat_centr_overlap$feature,
      start = check_genefeat_centr_overlap$start,
      end = check_genefeat_centr_overlap$end,
      score = check_genefeat_centr_overlap$score,
      strand = check_genefeat_centr_overlap$strand,
      frame = check_genefeat_centr_overlap$frame,
      attached = check_genefeat_centr_overlap$attribute,
      unique_ID = check_genefeat_centr_overlap$unique_ID,
      CentrStart = new_split_coord_for_overlap_genes_with_centrom,
      ChrEND = check_genefeat_centr_overlap$ChrEND
    )
    # 5.3 Write out new split coordinates for genes overlapping centromeres
    # write.table(new_split_coord_for_overlap_genes_with_centrom,
    #    file = paste0(out_path, "test_overlap.bed"),
    #    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    # 5.4 Prepare final BED with only Chr, CentrStart, ChrEND columns
    BEDnew_split_coord_for_overlap_genes_with_centrom = data.frame(
      Chr = new_split_coord_for_overlap_genes_with_centrom$Chr,
      CentrStart = new_split_coord_for_overlap_genes_with_centrom$CentrStart,
      ChrEND = new_split_coord_for_overlap_genes_with_centrom$ChrEND
    ) %>%
      dplyr::distinct()

    # --------------------------------------------------------------------------
    # 6 Combine final BED coordinates
    # --------------------------------------------------------------------------
    # 6.1 Get rows for chromosomes without overlapping genes
    df_exact_centromere_split = read_cent[
      read_cent$Chr %in%
        BEDnew_split_coord_for_overlap_genes_with_centrom$Chr ==
        FALSE,
    ]
    # 6.2 Combine
    BED_df_combined = rbind(
      df_exact_centromere_split,
      BEDnew_split_coord_for_overlap_genes_with_centrom
    )
    # 6.3 Sort as original bed file
    BED_df_combined_sorted_as_original = BED_df_combined[
      match(read_cent$Chr, BED_df_combined$Chr),
    ]

    # for GTF
    write.table(
      BED_df_combined_sorted_as_original,
      file = paste0(out_path, "RECAL1_FINAL_SPLIT_coord_for_GTF.bed"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
    # for FASTA
    write.table(
      prep_centromere_BED(
        BED_df_combined_sorted_as_original,
        contig_prefix,
        mt_name,
        pt_name
      ),
      file = paste0(out_path, "RECAL1_FINAL_SPLIT_coord_for_FASTA.bed"),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  } else {
    #5.5 Write out BED with only Chr, CentrStart, ChrEND columns
    # write.table(BEDnew_split_coord_for_overlap_genes_with_centrom,
    #   file = paste0(out_path, "test.bed"),
    #   quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    message(
      " \n",
      "No overlapping gene features detected use the original centromere coordinates BED file.\n",
      "Use:\n",
      centromere_path,
      " as the input to modify GTF file.\n",
      "and:\n",
      out_path,
      "IN0_original_centromeres_for_FASTA.bed",
      " for splitting FASTA.",
    )
  }
}

#########################################################################################################################################
# how to run:
#prepare_centromere_split = function(gtf_path, centromere_path, out_path, contig_prefix, mt_name, pt_name)
# prepare_centromere_split(
#   gtf_path = "/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/A0_toy_all_scenarios_2chr.gtf,   #A0_toy_all_scenarios.gtf",     ###/A0_toy_all_scenarios.gtf",
#   centromere_path = "/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/A0_toy_centromeres_for_gtf.bed",
#   out_path = "/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/",
#    "CAJHDD",
#    "Mt",
#    "Pt"
# )
