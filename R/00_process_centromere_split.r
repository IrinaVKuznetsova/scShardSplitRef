#' Check for gene feature overlaps with centromere coordinates and prepare new split coordinates if needed
#'
#' Converts a BED file to a format suitable for splitting GTF files
#'
#' @param df A data frame containing centromere coordinates with columns "Chr",
#'  "CentrStart", and "ChrEND".
#' @param contig_prefix A string prefix to identify contigs in the chromosome
#'  names.
#' @param mt_name A string representing the name of the mitochondrial chromosome
#'  (*e.g.*, "Mt").
#' @param pt_name A string representing the name of the plastid chromosome
#'  (*e.g.*, "Pt").
#'
#' @importFrom utils read.table write.table
#'
#' @examples
#' prepare_centromere_BED(df,
#'                        contig_prefix = "CAJHDD",
#'                        mt_name = "Mt",
#'                        pt_name = "Pt")
#'
#' @author Irina Kuznetsova, \email{irina.Kuznetsova@@curtin.edu.au}
#' @dev
#'
################################################################################
# FUNCTION-1 checks for the rare case where a gene feature lies on the
# centromere coordinate and recalculates the split coordinate.
# It takes the end coordinate of any gene that overlaps the centromere and
# shifts it 10 bp downstream to establish a new split coordinate.
#
# FUNCTION-0 prepares centromere coordinates in the correct format for FASTA
# splitting.
################################################################################

prepare_centromere_BED = function(df, contig_prefix, mt_name, pt_name) {
  contig <- grepl(sprintf("^%s", contig_prefix), df$Chr)
  organelle <- df$Chr %in% c(mt_name, pt_name)
  to_keep <- contig | organelle

  split_rows <- (!to_keep) & (df$CentrStart > 0L)
  unchanged_rows <- to_keep | (!to_keep & df$CentrStart <= 0L)

  df_unchanged <- df[unchanged_rows, , drop = FALSE]

  if (any(split_rows)) {
    chr_split <- df$Chr[split_rows]
    start_split <- df$CentrStart[split_rows]
    end_split <- df$ChrEND[split_rows]

    df_split <- data.frame(
      Chr = rep(chr_split, each = 2L),
      CentrStart = as.vector(rbind(0L, start_split)),
      ChrEND = as.vector(rbind(start_split, end_split))
    )

    df_result <- rbind(df_unchanged, df_split)
    rownames(df_result) <- NULL
  } else {
    df_result <- df_unchanged
  }

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
prepare_centromere_split <- function(
  gtf_path,
  centromere_path,
  out_path,
  contig_prefix,
  mt_name,
  pt_name
) {
  # --------------------------------------------------
  # 1. Read GTF file, add column names
  # --------------------------------------------------
  gtf_file <- read.table(file = gtf_path, header = FALSE, sep = "\t")
  colnames(gtf_file) <- c(
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
  gtf_file$unique_ID <- sprintf("ID%d", seq_len(nrow(gtf_file)))

  # --------------------------------------------------
  # 2. Read centromere positions file (BED format)
  # --------------------------------------------------
  read_cent <- read.table(file = centromere_path, sep = "\t")
  colnames(read_cent) <- c("Chr", "CentrStart", "ChrEND")

  # Convert original BED format ready for FASTA processing
  orig_For_FASTA <- prepare_centromere_BED(
    read_cent,
    contig_prefix,
    mt_name,
    pt_name
  )
  write.table(
    orig_For_FASTA,
    file = sprintf("%sIN0_original_centromeres_for_FASTA.bed", out_path),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  # ---------------------------------------------------------------------------
  # 3. Combine the two data frames
  # --------------------------------------------------------------------------
  df_merged_chr <- merge(gtf_file, read_cent, by = "Chr")

  # --------------------------------------------------------------------------
  # 4. Check for overlapping gene features (vectorised)
  # --------------------------------------------------------------------------
  check_genefeat_centr_overlap <-
    df_merged_chr[
      df_merged_chr[["start"]] < df_merged_chr[["CentrStart"]] &
        df_merged_chr[["end"]] > df_merged_chr[["CentrStart"]],
    ]

  # Prepare empty object to store new split coordinates
  BEDnew_split_coord_for_overlap_genes_with_centrom <- data.frame()

  # --------------------------------------------------------------------------
  # 5. Handle overlaps
  # --------------------------------------------------------------------------
  if (nrow(check_genefeat_centr_overlap) > 0) {
    cli::cli_alert_warning(
      "{.strong Overlapping gene features with centromere coordinates detected!}
      Re-calculating split coordinates for these chromosomes.

      Use:
        {.file {sprintf('%sRECAL1_FINAL_SPLIT_coord_for_GTF.bed', out_path)}}
        as the input to modify the GTF file; and

        {.file {sprintf('%sRECAL1_FINAL_SPLIT_coord_for_FASTA.bed', out_path)}}
        as the input for splitting FASTA."
    )

    # 5.1 Shift centromere split coordinate 10 bp downstream of the gene feature end coordinate
    new_split_coord <- check_genefeat_centr_overlap$end + 10L

    # 5.2 Create new data frame with updated split coordinates
    new_overlap_df <- data.frame(
      Chr = check_genefeat_centr_overlap$Chr,
      source = check_genefeat_centr_overlap$source,
      feature = check_genefeat_centr_overlap$feature,
      start = check_genefeat_centr_overlap$start,
      end = check_genefeat_centr_overlap$end,
      score = check_genefeat_centr_overlap$score,
      strand = check_genefeat_centr_overlap$strand,
      frame = check_genefeat_centr_overlap$frame,
      attribute = check_genefeat_centr_overlap$attribute,
      unique_ID = check_genefeat_centr_overlap$unique_ID,
      CentrStart = new_split_coord,
      ChrEND = check_genefeat_centr_overlap$ChrEND
    )

    # 5.4 Prepare the final BED with only Chr, CentrStart, ChrEND columns and drop duplicates
    tmp_bed <- unique(new_overlap_df[, c("Chr", "CentrStart", "ChrEND")])

    # --------------------------------------------------------------------------
    # 6. Combine final BED coordinates
    # --------------------------------------------------------------------------
    # 6.1 Get rows for chromosomes without overlapping genes
    non_overlap_chr <- !(read_cent$Chr %in% tmp_bed$Chr)
    df_exact_centromere_split <- read_cent[non_overlap_chr, ]

    # 6.2 Combine
    BED_df_combined <- rbind(df_exact_centromere_split, tmp_bed)

    # 6.3 Sort as original BED file order
    order_idx <- match(read_cent$Chr, BED_df_combined$Chr)
    BED_df_combined_sorted_as_original <- BED_df_combined[order_idx, ]

    gtf_bed_outfile <- sprintf(
      "%sRECAL1_FINAL_SPLIT_coord_for_GTF.bed",
      out_path
    )
    fasta_bed_outfile <- sprintf(
      "%sRECAL1_FINAL_SPLIT_coord_for_FASTA.bed",
      out_path
    )
    # for GTF
    write.table(
      BED_df_combined_sorted_as_original,
      file = gtf_bed_outfile,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
    # for FASTA
    write.table(
      prepare_centromere_BED(
        BED_df_combined_sorted_as_original,
        contig_prefix,
        mt_name,
        pt_name
      ),
      file = fasta_bed_outfile,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  } else {
    cli::cli_alert_success(
      "No overlapping gene features detected.
      Use:
        {.file {centromere_path}}
      as the input to modify GTF file.
      and:
        {.file {sprintf('%sIN0_original_centromeres_for_FASTA.bed', out_path)}}
      for splitting FASTA."
    )
  }
}
