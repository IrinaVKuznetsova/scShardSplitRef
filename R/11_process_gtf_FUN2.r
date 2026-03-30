#' Process the GTF file
#'
#' @param gtf_path A string representing the file path to the GTF file.
#' @param centromere_path A string representing the file path to the centromere
#'   coordinates BED file.
#' @param unchar_region TODO: describe this param.
#' @param mito TODO: describe this param.
#' @param chloro TODO: describe this param.
#'
#' @examples
#' #TODO: Add examples of use for this fn.
#'
#' @inheritParams prepare_centromere_split
#' @autoglobal

process_gtf <- function(
  gtf_path,
  centromere_path,
  unchar_region,
  mito,
  chloro
) {
  # Validation and reading
  cli::cli_alert_info("Validating inputs and reading files...")
  validated <- validate_inputs(gtf_path, centromere_path)
  gtf_file <- validated$gtf
  read_cent <- validated$bed

  # Add unique ID for ordering
  gtf_file$unique_ID <- sprintf("ID%i", seq_len(nrow(gtf_file)))

  # Merge GTF/BED by chromosome
  df_merged_chr <- merge(gtf_file, read_cent, by = "Chr", all.x = TRUE)

  # Build regex for special chromosomes
  contig_regex <- sprintf("^(%s|%s|%s)$", unchar_region, mito, chloro)

  cli::cli_alert_info("Splitting contigs and chromosomes...")

  # Identify special contigs (e.g., CAJHD, Mt, Pt)
  is_special <- grepl(contig_regex, df_merged_chr$Chr)

  contigs_rows <- df_merged_chr[is_special, , drop = FALSE]
  contigs_rows$NEW <- sprintf(
    "%s-%s-%s",
    contigs_rows$Chr,
    contigs_rows$CentrStart,
    contigs_rows$ChrEND
  )

  chr_df <- df_merged_chr[!is_special, , drop = FALSE]

  # Split by gene position relative to centromere
  # 1. Below centromere
  idx_less_than <- with(chr_df, start < CentrStart & end <= CentrStart)
  chr_less_than <- chr_df[idx_less_than, , drop = FALSE]
  chr_less_than$NEW <- sprintf(
    "%s-%d-%s",
    chr_less_than$Chr,
    0L,
    chr_less_than$CentrStart
  )

  # 2. Above centromere
  idx_greater_than <- with(chr_df, start >= CentrStart & end > CentrStart)
  chr_greater_than <- chr_df[idx_greater_than, , drop = FALSE]
  chr_greater_than$NEW <- sprintf(
    "%s-%s-%s",
    chr_greater_than$Chr,
    chr_greater_than$CentrStart,
    chr_greater_than$ChrEND
  )

  # 3. Overlaps centromere
  idx_overlap <- with(chr_df, start < CentrStart & end > CentrStart)
  chr_overlap <- chr_df[idx_overlap, , drop = FALSE]
  chr_overlap$NEW <- sprintf(
    "%s-%s-%s",
    chr_overlap$Chr,
    chr_overlap$CentrStart,
    chr_overlap$ChrEND
  )
  if (nrow(chr_overlap) > 0L) {
    cli::cli_warn(c(
      "Some gene features overlap the centromere boundary.",
      "i" = "{nrow(chr_overlap)} detected. This can indicate BED/GTF
      inconsistency.",
      "i" = "See 'prepare_centromere_split.r' for guidance."
    ))
  }

  # Recalculate coordinates for genes above centromere
  cli::cli_alert_info(
    "Recalculating coordinates for features above centromere..."
  )
  feature_length <- chr_greater_than$end - chr_greater_than$start
  dist_chr_end_feat_end <- chr_greater_than$ChrEND - chr_greater_than$end
  dist_centr_chr_end <- chr_greater_than$ChrEND - chr_greater_than$CentrStart

  chr_greater_than$start <- dist_centr_chr_end -
    dist_chr_end_feat_end -
    feature_length
  chr_greater_than$end <- dist_centr_chr_end - dist_chr_end_feat_end

  # Combine all subsets and restore input order by unique ID
  cli::cli_alert_info("Combining all feature rows...")
  combined <- rbind(chr_greater_than, chr_less_than, chr_overlap, contigs_rows)
  combined <- combined[
    match(gtf_file$unique_ID, combined$unique_ID),
    ,
    drop = FALSE
  ]

  # Form GTF output columns
  out <- combined[, c(
    "NEW",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )]

  cli::cli_alert_success("Processing complete. Returning formatted GTF.")

  return(out)
}
