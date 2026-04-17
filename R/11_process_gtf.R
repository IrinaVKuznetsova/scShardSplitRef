#' Process the GTF file
#'
#' @param gtf_path A string representing the file path to the GTF file.
#' @param centromere_path A string representing the file path to the centromere
#'   coordinates BED file.
#' @param unchar_region A string specifying the pattern for uncharacterized
#'   region chromosome names.
#' @param mito A string specifying the mitochondrial chromosome name.
#' @param chloro A string specifying the chloroplast chromosome name.
#'
#' @examples
#'
#' gtf <- system.file("extdata",
#'          "A3_toy_all_scenarios_2chr.gtf",
#'          package = "scShardSplitRef",
#'          mustWork = TRUE
#'          )
#' bed <- system.file("extdata",
#'                    "IN0_toy_centromeres_for_gtf.bed",
#'                    package = "scShardSplitRef",
#'                    mustWork = TRUE)
#' result <- process_gtf(
#'    gtf_path = gtf,
#'    centromere_path = bed,
#'    unchar_region = "CAJHDD.*",
#'    mito = "Mt",
#'    chloro = "Pt"
#'    )
#'
#' @returns A data frame with formatted GTF columns.
#' @export
#' @autoglobal

process_gtf <- function(
  gtf_path,
  centromere_path,
  unchar_region,
  mito,
  chloro
) {
  # Validation and reading
  cli::cli_alert("Validating inputs and reading files...")
  validated <- validate_inputs(gtf_path, centromere_path)
  gtf_file <- validated$gtf
  bed_file <- validated$bed

  # Add unique ID for ordering
  gtf_file$unique_ID <- sprintf("ID%i", seq_len(nrow(gtf_file)))

  # Standardize BED column names for internal use
  bed_file$CentrStart <- bed_file$RegionStart
  bed_file$ChrEND <- bed_file$RegionEnd

  # Merge GTF/BED by chromosome
  df_merged_chr <- merge(
    gtf_file,
    bed_file[, c("Chr", "CentrStart", "ChrEND")],
    by = "Chr",
    all.x = TRUE
  )

  # Build regex for special chromosomes
  contig_regex <- sprintf("^(%s|%s|%s)$", unchar_region, mito, chloro)

  cli::cli_alert("Splitting contigs and chromosomes...")

  # Identify special contigs (e.g., CAJHD, Mt, Pt)
  is_special <- grepl(contig_regex, df_merged_chr$Chr)
  contigs_rows <- .assign_coords(
    df_merged_chr[is_special, , drop = FALSE],
    df_merged_chr[is_special, "Chr"],
    df_merged_chr[is_special, "CentrStart"],
    df_merged_chr[is_special, "ChrEND"]
  )

  chr_df <- df_merged_chr[!is_special, , drop = FALSE]

  # Remove rows with NA centromere info from chr_df
  chr_df <- chr_df[!is.na(chr_df$CentrStart), , drop = FALSE]

  # Define indices for chromosome position relative to centromere
  idx_less_than <- with(chr_df, start < CentrStart & end <= CentrStart)
  idx_greater_than <- with(chr_df, start >= CentrStart & end > CentrStart)
  idx_overlap <- with(chr_df, start < CentrStart & end > CentrStart)

  chr_less_than <- .assign_coords(
    chr_df[idx_less_than, , drop = FALSE],
    chr_df[idx_less_than, "Chr"],
    0L,
    chr_df[idx_less_than, "CentrStart"]
  )

  chr_greater_than <- .assign_coords(
    chr_df[idx_greater_than, , drop = FALSE],
    chr_df[idx_greater_than, "Chr"],
    chr_df[idx_greater_than, "CentrStart"],
    chr_df[idx_greater_than, "ChrEND"]
  )

  chr_overlap <- .assign_coords(
    chr_df[idx_overlap, , drop = FALSE],
    chr_df[idx_overlap, "Chr"],
    chr_df[idx_overlap, "CentrStart"],
    chr_df[idx_overlap, "ChrEND"]
  )

  if (nrow(chr_overlap) > 0L) {
    cli::cli_warn(c(
      "Some gene features overlap the centromere boundary.",
      x = "{nrow(chr_overlap)} detected. This can indicate BED/GTF inconsistency.",
      i = "See 'prepare_centromere_split.r' for guidance."
    ))
  }

  # Recalculate coordinates for genes above centromere
  cli::cli_alert(
    "Recalculating coordinates for features above centromere..."
  )
  chr_greater_than <- .recalculate_above_centromere_coords(chr_greater_than)

  # Combine all subsets and restore input order by unique ID
  cli::cli_alert("Combining all feature rows...")
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

  cli::cli_alert("Processing complete. Returning formatted GTF.")

  return(out)
}

#' Assign NEW coordinate column to features
#'
#' @param df A data frame to receive the NEW column.
#' @param chr Character vector for Chr component.
#' @param start Numeric vector for start component.
#' @param end Numeric vector for end component.
#'
#' @returns The input data frame with NEW column added.
#' @dev

.assign_coords <- function(df, chr, start, end) {
  if (nrow(df) > 0L) {
    df$NEW <- sprintf("%s-%s-%s", chr, start, end)
  }
  df
}

#' Recalculate start and end coordinates for features above centromere
#'
#' @param df A data frame of features above the centromere with start, end,
#'   ChrEND, CentrStart columns.
#'
#' @returns The input data frame with recalculated start and end coordinates.
#' @dev

.recalculate_above_centromere_coords <- function(df) {
  feature_length <- df$end - df$start
  dist_chr_end_feat_end <- df$ChrEND - df$end
  dist_centr_chr_end <- df$ChrEND - df$CentrStart

  df$start <- dist_centr_chr_end - dist_chr_end_feat_end - feature_length
  df$end <- dist_centr_chr_end - dist_chr_end_feat_end

  df
}
