# Oct 2026
# by IK


# Genomes tested:
#   * Barley (Hordeum vulgare, Morex3)        - starts at line: 20-120
#   * Bread wheat (Triticum aestivum, Renan)  - starts at line: 126-221
#   * Faba bean (Vicia faba, Hedin2)          - starts at line: 226-325
#
# Note, Wheat and faba bean were built on the HPC scratch filesystem
#
# Validation workflow for scShardSplitRef
# This file documents both R and Bash commands.
# Step 1: Run in RStudio (R).
# Step 2: Run in a Bash terminal (bedtools).
# Step 3: Run in a Bash terminal (Cell Ranger ARC mkref).
# Do not execute the entire file as an R script.



############################################################
############################################################
# I. Barley (Morex v3) (24 July 2026)
############################################################
############################################################

#################################
# 1. Run in RStudio (R)
#################################
# 1.1 Use scShardSplitRef package
# PATHs
gtf_path="/path/to/input/01_Hordeum_vulgare.Mt_Pt_v49_Barley_v62_COMBINED.gtf"
bed_path="/path/to/input/02_GTF_FORMATTED_MorexV3_centromere_positions_MANUALLY_ALL.bed"

output_file_name_dir="/path/to/output/0_Barley/"


# 1.2 Determine split regions
get_split_reg <- determine_split_regions(
  bed = bed_path,
  gtf = gtf_path,
  output_bed = file.path(output_file_name_dir, "1_barley_split_regions.bed"),
  limit = 2L^29L,
  shift_by = 2000,
  clearance = 900
)  # the output is file


# 1.3 Modify GTF file
process_gtf(
  split_regions_bed = get_split_reg,
  gtf = gtf_path,
  genome_name = "subset_barley",
  genome_version = "3",
  out_path = output_file_name_dir
)


#####################################
# 2. Run in Bash terminal (bedtools)
####################################
# 2.1 Prep FASTA
# conda activate bio  # using conda environment 'bio' with required tools installed

IN_FA="/path/to/input/01_Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly_COMBINED_Mt_Pt_v49_Barley_v62.fa"
BED_FIN_OUT="/path/to/output/0_Barley/1_barley_split_regions.bed"
FASTA_FIN_OUT="/path/to/output/0_Barley/"

bedtools getfasta -fi ${IN_FA} -bed ${BED_FIN_OUT} -fo ${FASTA_FIN_OUT}/morex3_SPLIT.fasta

# Output, barley
# >1H:0-206486643
# >1H:206486643-516505932
# >2H:0-301293086
# >2H:301293086-665585731
# >3H:0-267852507
# >3H:267852507-621516506
# >4H:0-276149121
# >4H:276149121-610333535


# 2.2 Convert FASTA headers from `Chr:RegionStart-RegionEnd` to `Chr-RegionStart-RegionEnd`
sed '/^>/s/:/-/g' ${FASTA_FIN_OUT}/morex3_SPLIT.fasta >${FASTA_FIN_OUT}/morex3_split_FINAL.fasta

# 2.3 Validate
grep ">" ${FASTA_FIN_OUT}/morex3_split_FINAL.fasta

# Output
# >1H-0-206486643
# >1H-206486643-516505932
# >2H-0-301293086
# >2H-301293086-665585731
# >3H-0-267852507
# >3H-267852507-621516506
# >4H-0-276149121


#################################################
# 3. Run in Bash terminal (Cell Ranger ARC)
##############################################
# 3.1 Create CONFIG file as per Cell Ranger ARC instructions
# CONFIG file:  barley_mkref.config
# BASH SCRIPT:  0_barley_CellRangerARC.sh
INCONFIG=/path/to/output/0_Barley/
  
  cat > "${INCONFIG}/barley_mkref.config" <<'EOF'
{
  organism: "barley"
  genome: ["MOREX3_v3_Oct26"]
  input_fasta: ["/path/to/output/0_Barley/morex3_split_FINAL.fasta"]
  input_gtf: ["/path/to/output/0_Barley/B1_FINAL_MODIFIED_GTF_subset_barley_3.gtf"]
}
EOF

# 3.2 Run Cell Ranger ARC mkref
cellranger-arc mkref \
--config="${INCONFIG}/barley_mkref.config" \
--nthreads=24 \
> "${INCONFIG}/barley_mkref_new.log" 2>&1

tail -n 30 "${INCONFIG}/barley_mkref_new.log"





############################################################
############################################################
# II. Bread wheat (Triticum_aestivum_renan) (24 July 2026)
############################################################
############################################################

#################################
# 1. Run in RStudio (R)
#################################
# 1.1 Use scShardSplitRef package
# PATHs
gtf_path="/path/to/input/Triticum_aestivum_renan.Triticum_aestivum_Renan_v2.1.62.gtf"
bed_path="/path/to/input/wheat_chromosomes.bed"  # as no centromere info availble,  use chr length

output_file_name_dir="/path/to/output/1_Wheat/"


# 1.2 Determine split regions
get_split_reg <- determine_split_regions(
  bed = bed_path,
  gtf = gtf_path,
  output_bed = file.path(output_file_name_dir, "1_wheat_split_regions.bed"),
  limit = 2L^29L,
  shift_by = 2000,
  clearance = 900
)  # the output is file


# 1.3 Modify GTF file
process_gtf(
  split_regions_bed = get_split_reg,
  gtf = gtf_path,
  genome_name = "wheat",
  genome_version = "3",
  out_path = output_file_name_dir
)


#####################################
# 2.Run in Bash terminal (bedtools)
#####################################
# 2.1 Prep FASTA
# conda activate bio  # using conda environment 'bio' with required tools installed

IN_FA="/path/to/input/Triticum_aestivum_renan.Triticum_aestivum_Renan_v2.1.dna.toplevel.fa"
BED_FIN_OUT="/path/to/output/1_Wheat/1_wheat_split_regions.bed"
FASTA_FIN_OUT="/path/to/output/1_Wheat/"

bedtools getfasta -fi ${IN_FA} -bed ${BED_FIN_OUT} -fo ${FASTA_FIN_OUT}/wheat_SPLIT.fasta

# 2.2 Validate
grep ">" ${FASTA_FIN_OUT}/wheat_SPLIT.fasta

# Output, wheat
# >1A:0-296965174
# >1A:296965174-593930347
# >1B:0-351387832
# >1B:351387832-702775664
# >1D:0-494594617


# 2.3 Convert FASTA headers from `Chr:RegionStart-RegionEnd` to `Chr-RegionStart-RegionEnd`
sed '/^>/s/:/-/g' ${FASTA_FIN_OUT}/wheat_SPLIT.fasta >${FASTA_FIN_OUT}/wheat_split_FINAL.fasta

# 2.4 Validate
grep ">" ${FASTA_FIN_OUT}/wheat_split_FINAL.fasta

# Output
# >1A-0-296965174
# >1A-296965174-593930347
# >1B-0-351387832
# >1B-351387832-702775664
# >1D-0-494594617


###################################################
# 3. Run in Bash terminal (Cell Ranger ARC)
################################################
# 3.1 Create CONFIG file as per Cell Ranger ARC instructions

INCONFIG=/path/to/output/1_Wheat/
  
  cat > "${INCONFIG}/1_wheat_mkref.config" <<'EOF'
{
  organism: "wheat"
  genome: ["Renan_v3_Oct26"]
  input_fasta: ["/path/to/output/1_Wheat/wheat_split_FINAL.fasta"]
  input_gtf: ["/path/to/output/1_Wheat/B1_FINAL_MODIFIED_GTF_wheat_3.gtf"]
}
EOF

# 3.2 Run Cell Ranger ARC mkref
cellranger-arc mkref \
--config="${INCONFIG}/1_wheat_mkref.config" \
--nthreads=24 \
> "${INCONFIG}/1_wheat_mkref_new.log" 2>&1




############################################################
############################################################
# III. Faba bean (Hedin2_genome_v1.63) (24 July 2026)
############################################################
############################################################

#################################
# 1. Run in RStudio (R)
#################################
# 1.1 Use scShardSplitRef package
# PATHs
gtf_path="/path/to/input/Vicia_faba.Hedin2_genome_v1.63.gtf"
bed_path="/path/to/input/faba_bean_chromosomes.bed"

output_file_name_dir="/path/to/output/2_FabaBean/"


# 1.2 Determine split regions
get_split_reg <- determine_split_regions(
  bed = bed_path,
  gtf = gtf_path,
  output_bed = file.path(output_file_name_dir, "1_fababean_split_regions.bed"),
  limit = 2L^29L,
  shift_by = 2000,
  clearance = 900
)  # the output is file


# 1.3 Modify GTF file
process_gtf(
  split_regions_bed = get_split_reg,
  gtf = gtf_path,
  genome_name = "Hedin2",
  genome_version = "3",
  out_path = output_file_name_dir
)


######################################
# 2. Run in Bash terminal (bedtools)
######################################
# 2.1 Prep FASTA
# conda activate bio  # using conda environment 'bio' with required tools installed

IN_FA="/path/to/input/vf_combined.fa"
BED_FIN_OUT="/path/to/output/2_FabaBean/1_fababean_split_regions.bed"
FASTA_FIN_OUT="/path/to/output/2_FabaBean/"

bedtools getfasta -fi ${IN_FA} -bed ${BED_FIN_OUT} -fo ${FASTA_FIN_OUT}/fababean_SPLIT.fasta

# 2.2 Validate
grep ">" ${FASTA_FIN_OUT}/fababean_SPLIT.fasta

# Output, faba bean
# >OX451735.1:0-524842364
# >OX451735.1:524842364-1049684729
# >OX451735.1:1049684729-1574527093
# >OX451736.1:0-451311207
# >OX451736.1:451311207-902622414
# >OX451736.1:902622414-1353933622
# >OX451736.1:1353933622-1805244829
# >OX451737.1:0-429192404
# >OX451737.1:429192404-858384808
# >OX451737.1:858384808-1287577211
# >OX451737.1:1287577211-1716769615
# >OX451738.1:0-409453994


# 2.3 Convert FASTA headers from `Chr:RegionStart-RegionEnd` to `Chr-RegionStart-RegionEnd`
sed '/^>/s/:/-/g' ${FASTA_FIN_OUT}/fababean_SPLIT.fasta >${FASTA_FIN_OUT}/fababean_split_FINAL.fasta

# 2.4 Validate
grep ">" ${FASTA_FIN_OUT}/fababean_split_FINAL.fasta

# Output
# >OX451735.1-524842364-1049684729
# >OX451735.1-1049684729-1574527093
# >OX451736.1-0-451311207


###################################################
# 3. Run in Bash terminal (Cell Ranger ARC)
################################################
# 3.1 Create CONFIG file as per Cell Ranger ARC instructions
INCONFIG=/path/to/output/2_FabaBean/
  
  cat > "${INCONFIG}/2_fababean_CR_ARC.config" <<'EOF'
{
  organism: "fababean"
  genome: ["Hedin2_Oct26"]
  input_fasta: ["/path/to/output/2_FabaBean/fababean_split_FINAL.fasta"]
  input_gtf: ["/path/to/output/2_FabaBean/B1_FINAL_MODIFIED_GTF_Hedin2_3.gtf"]
}
EOF

# 3.2 Run Cell Ranger ARC mkref
cellranger-arc mkref \
--config="${INCONFIG}/2_fababean_CR_ARC.config" \
--nthreads=24 \
> "${INCONFIG}/1_fababean_mkref_new.log" 2>&1
