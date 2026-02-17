#' @import dplyr
#' @import stringr
#' @importFrom utils read.table write.table













# 16 Feb 2026
# IK
# License:  GPL-3.0 



# FUNCTION-2
# This function takes GTF and BED (centromere) files as input and modifies the 1st, 5th, and 6th columns 
# according to the position of each gene feature on the chromosome and its relation to the centromere coordinates.  


# mamba activate 00_R_env
#library(tidyverse)
process_gtf = function(gtf_path, centromere_path, unchar_region, mito, chloro) {
#gtf_path="/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_input/A0_toy_all_scenarios_2chr.gtf"
#centromere_path="/data/02_PROJECTS/sc_ShardSplitRef/01_Working_R_scripts/Toy_validate_GTF_script/toy_output/RECAL1_FINAL_SPLIT_coord_for_GTF.bed"
#unchar_region="CAJHD"
#mito="Mt"
#chloro="Pt"
    contig_regex = paste0( "(", unchar_region, "|", mito, "|", chloro, ")" )

# --- FILE CHECKS ---------------------------------------------------------
    validated = validate_inputs(gtf_path, centromere_path)
    gtf_file = validated$gtf
    read_cent = validated$bed

    # Add unique ID
    gtf_file$unique_ID = paste0("ID", seq_len(nrow(gtf_file)))

# --- END (FILE CHECKS) END -------------------------------------------------------------

    # ---------------------------------------------------------------------------
    # 3 Combine two data frames (GTF file and centromere position by 'Chr')
    df_merged_chr = merge(gtf_file, read_cent, by = "Chr")

# NOTE, barley geneome has 1-7 chromosome, but other might be larger and also have different contig suffix and Mt, Pt annotations
# thus, this can go to function parameters !!!!! will need to add this later (4.1 and 4.2)

    # ---------------------------------------------------------------------------
    # 4.1 Subset df to CAJHDD, Mt, Pt # str_detect(Chr, '^CAJHDD|^Mt$|^Pt$')  ### ensures matches to CAJHDD contigs only; Mt exactly; Pt exactly
    # ***this subset will be attached to the dataframe with re-annotated 1st column with coordinates deppending on the centromere, to keep all values unique 
    contigs_rows = df_merged_chr %>%
        filter(str_detect(Chr, contig_regex)) %>%
        mutate(NEW= paste( Chr, CentrStart, ChrEND, sep="-"))
    # -----------------------------------------------------------------------------
    # 4.2 Subset df to everything except CAJHDD, Mt, Pt   ////"CAJHDD01"/ chrUn
    chr_df = df_merged_chr %>%
        filter(!str_detect(Chr, contig_regex)) 
    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    
    # TEST where a gene feature cooridnate lays BELOW the centromere position or ABOVE
    # 4.3  Subset df[1H-7H] with condition that   gene features 'START-END' coordinate is less than "CENTROMERE" coordinate | less: 1H-0-206486643
    # this accounts for both "start" and "end" gene features being within 0 - centromere position 
    chr_less_than = chr_df %>%
        #filter(end < CentrStart ) %>%   ## check
        filter(start < CentrStart & end <= CentrStart) %>%
        mutate(NEW= paste( Chr, 0, CentrStart, sep="-"))
# (start < CentrStart & end <= CentrStart) - covers following cases:
# BED: 1H   100   200
# GTF:
# 1H  source  gene    0   20  .   +   .   gene_id "gene1"; gene_name "gene1";
# 1H  source  transcript  0   20  .   +   .   gene_id "gene1"; transcript_id "gene1.t1";
# 1H  source  exon    0   10  .   +   .   gene_id "gene1"; transcript_id "gene1.t1"; exon_number "1";
# 1H  source  exon    11  20  .   +   .   gene_id "gene1"; transcript_id "gene1.t1"; exon_number "2";
# 1H  source  gene    20  40  .   -   .   gene_id "gene2"; gene_name "gene2";
# 1H  source  transcript  20  40  .   -   .   gene_id "gene2"; transcript_id "gene2.t1";
# 1H  source  exon    20  28  .   -   .   gene_id "gene2"; transcript_id "gene2.t1"; exon_number "1";
# 1H  source  exon    29  40  .   -   .   gene_id "gene2"; transcript_id "gene2.t1"; exon_number "2";
# 1H  source  gene    80  100 .   +   .   gene_id "gene3"; gene_name "gene3"; note "ends_at_centromere";
# 1H  source  transcript  80  100 .   +   .   gene_id "gene3"; transcript_id "gene3.t1";
# 1H  source  exon    80  90  .   +   .   gene_id "gene3"; transcript_id "gene3.t1"; exon_number "1";
# 1H  source  exon    91  100 .   +   .   gene_id "gene3"; transcript_id "gene3.t1"; exon_number "2";


    # ------------------------------------------------------------------------------------------------------------------------------------------
    # 4.4 Subset df[1H-7H] with condition that gene feature 'END' coordinate is greater than "CENTROMERE" greater:  1H-206486643-516505932
    chr_greater_than = chr_df %>%
        filter(start >= CentrStart & end > CentrStart) %>%  
        ##filter( end > CentrStart ) %>%
        mutate(NEW= paste( Chr, CentrStart, ChrEND, sep="-"))
# (start >= CentrStart & end > CentrStart) - covers following cases:
# Gene features for this scenario needs to be modified
# BED: 1H   100   200
# GTF:
# 1H  source  gene    100 120 .   +   .   gene_id "gene4"; gene_name "gene4"; note "starts_at_centromere";
# 1H  source  transcript  100 120 .   +   .   gene_id "gene4"; transcript_id "gene4.t1";
# 1H  source  exon    100 110 .   +   .   gene_id "gene4"; transcript_id "gene4.t1"; exon_number "1";
# 1H  source  exon    111 120 .   +   .   gene_id "gene4"; transcript_id "gene4.t1"; exon_number "2";
# 1H  source  gene    120    140    .   -   .   gene_id "gene5"; gene_name "gene5";
# 1H  source  transcript  120    140    .   -   .   gene_id "gene5"; transcript_id "gene5.t1";
# 1H  source  exon    120    130    .   -   .   gene_id "gene5"; transcript_id "gene5.t1"; exon_number "1";
# 1H  source  exon    131    140    .   -   .   gene_id "gene5"; transcript_id "gene5.t1"; exon_number "2";
# 1H  source  gene    180 200 .   +   .   gene_id "gene6"; gene_name "gene6"; note "ends_at_chromosome_end";
# 1H  source  transcript  180 200 .   +   .   gene_id "gene6"; transcript_id "gene6.t1";
# 1H  source  exon    180 190 .   +   .   gene_id "gene6"; transcript_id "gene6.t1"; exon_number "1";
# 1H  source  exon    191 200 .   +   .   gene_id "gene6"; transcript_id "gene6.t1"; exon_number "2";



# ------------------------------------------------------------------------------------------------------------------------------------------
    # 4.5 Handle overlapping features - this to double-check in case a wrong centromere position file is used
    # (features that cross the centromere boundary).  
    chr_overlap = chr_df %>%
        #filter(start < CentrStart & end > CentrStart) %>%
        filter(start < CentrStart & end > CentrStart) %>%
        mutate(NEW = paste(Chr, CentrStart, ChrEND, sep="-"))
    cat(
    "Number of rows:", nrow(chr_overlap), "\n",
    "This should be zero if you are using the correct centromere BED file.\n",
    "If the number of rows is greater than 0, make sure you have run the script prepare_centromere_split.r\n",
    "which checks for rare cases where a gene overlaps the centromere position\n",
    "and recalculates the centromere coordinates accordingly.\n"
    )


    # ------------------------------------------------------------------------------------------------------------------------------------------
    # 5 Recalculate start-end coordinates for features laying in CENTROMERE-ChrEND
    # 5.1 Calculate feature length [end-start]
    feature_length = chr_greater_than$end - chr_greater_than$start
    # Calculate distance from chromosome end to a feature-end:
    distan_ChrEnd_FeatEnd = chr_greater_than$ChrEND - chr_greater_than$end
    # Calculate distance from a feature-start to centromer:
    distan_Centr_FeatStart = chr_greater_than$start - chr_greater_than$CentrStart
    # Calculate distance from Centromere to Chr end
    distan_Centr_ChrEnd = chr_greater_than$ChrEND - chr_greater_than$CentrStart
    # NewFeatStartCoordint
    NewStartGene = distan_Centr_ChrEnd - distan_ChrEnd_FeatEnd - feature_length
    # NewFeatEndCoordint
    NewEndGene = distan_Centr_ChrEnd - distan_ChrEnd_FeatEnd 
    # 5.2 Append
    chr_greater_than$NewFeatStartCoordint = NewStartGene
    chr_greater_than$NewFeatEndCoordint = NewEndGene
    chr_greater_than$NewChrEnd = distan_Centr_ChrEnd
    chr_greater_than$NewChrStart = 0
    chr_greater_than$FeatureLength=feature_length
    chr_greater_than$distan_ChrEnd_FeatEnd = distan_ChrEnd_FeatEnd
    chr_greater_than$distan_Centr_FeatStart = distan_Centr_FeatStart
    # ------------------------------------------------------------------------------------------------------------------------------------------
    # 6 COMBINE all data into one
    # 6.1 Replace values in "start" and "end" with newly calculated coordinates:
    chr_greater_than$start = chr_greater_than$NewFeatStartCoordint
    chr_greater_than$end = chr_greater_than$NewFeatEndCoordint
    # 6.2 Remove duplicated columns
    chr_greater_than$NewFeatStartCoordint=NULL
    chr_greater_than$NewFeatEndCoordint=NULL
    chr_greater_than$FeatureLength= NULL
    chr_greater_than$NewChrEnd = NULL
    chr_greater_than$NewChrStart = NULL
    chr_greater_than$distan_ChrEnd_FeatEnd = NULL
    chr_greater_than$distan_Centr_FeatStart = NULL
    # 6.3 Put all obtained info into data frame
    FINAL_df = rbind(chr_greater_than, chr_less_than, chr_overlap, contigs_rows)   # added output from 4.5 
#    FINAL_df=rbind(chr_greater_than, chr_less_than, contigs_rows)
    df_ordered =  FINAL_df[ match(gtf_file$unique_ID, FINAL_df$unique_ID),  ] 
    # FIND in GTF format (9 columns)
    FINAL_df2= df_ordered[, c("NEW", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")]
    return(FINAL_df2)
}
