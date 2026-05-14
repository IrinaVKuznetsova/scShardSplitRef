# 2 April 2026
# IK


# Only a very small number of plant genomes contain any chromosome or contig larger than 536.8 Mb (2^29 bases).
#  Triticum aestivum (bread wheat)	IWGSC RefSeq v1.0/v2.1	Most chromosomes 600 -850 Mb
#  Hordeum vulgare (barley)	Morex V3/V4	Yes -several	2H, 3H, 5H exceed threshold
#  Avena sativa (oat)	Hexaploid; several >600 Mb
#  Secale cereale (rye)	Large 7-chromosome genome


#------------------------------------------------------------------------------------------------------
# 1.0 Triticum aestivum (bread wheat) - [15-17 Gigabases]
#   Genome assembly: Triticum_aestivum_Renan_v2.1
#   https://plants.ensembl.org/Triticum_aestivum_renan/Info/Index
#------------------------------------------------------------------------------------------------------
#cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/
mkdir 00_TriticumAestivumRenan
cd 00_TriticumAestivumRenan
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/triticum_aestivum_renan/dna/Triticum_aestivum_renan.Triticum_aestivum_Renan_v2.1.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gtf/triticum_aestivum_renan/Triticum_aestivum_renan.Triticum_aestivum_Renan_v2.1.62.gtf.gz	



#------------------------------------------------------------------------------------------------------
# 1.1 Secale cereale (rye) - [6.74 Gigabases]
#   Genome assembly: Rye_Lo7_2018_v1p1p1
#   https://plants.ensembl.org/Secale_cereale/Info/Index
#------------------------------------------------------------------------------------------------------
#cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/
mkdir 01_SecaleCereale
cd 01_SecaleCereale
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/secale_cereale/dna/Secale_cereale.Rye_Lo7_2018_v1p1p1.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gtf/secale_cereale/Secale_cereale.Rye_Lo7_2018_v1p1p1.62.gtf.gz	



#------------------------------------------------------------------------------------------------------
# 1.2 Avena_sativa (oat) - [~11 Gigabases]
#   Genome assembly: Oat_OT3098_v2
#   https://plants.ensembl.org/Avena_sativa_OT3098/Info/Index
#------------------------------------------------------------------------------------------------------
#cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/
mkdir 02_AvenaSativaOats
cd 02_AvenaSativaOats
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/avena_sativa_ot3098/dna/Avena_sativa_ot3098.Oat_OT3098_v2.dna.toplevel.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gtf/avena_sativa_ot3098/Avena_sativa_ot3098.Oat_OT3098_v2.62.gtf.gz



| Species              | Chr | Size (bp)     |
|----------------------|-----|---------------|
| **Exceeds Cell Ranger Limit** |     | 536,800,000 |
| **Oats**             | 2C  | 589,118,817   |
|                      | 3C  | 638,425,132   |
|                      | 4C  | 716,105,986   |
|                      | 5C  | 613,160,974   |
|                      | 6C  | 626,220,839   |
|                      | 7C  | 551,718,542   |
| **Secale cereale (rye)** | 1R  | 727,344,967   |
|                      | 2R  | 946,003,158   |
|                      | 3R  | 965,754,312   |
|                      | 4R  | 906,459,801   |
|                      | 5R  | 876,148,008   |
|                      | 6R  | 885,153,844   |
|                      | 7R  | 899,925,126   |
| **Wheat**            | 1A  | 593,930,347   |
|                      | 1B  | 702,775,664   |
|                      | 2A  | 792,837,209   |
|                      | 2B  | 812,232,696   |
|                      | 2D  | 661,835,603   |
|                      | 3A  | 750,337,041   |
|                      | 3B  | 854,463,248   |
|                      | 3D  | 623,248,023   |
|                      | 4A  | 749,950,614   |
|                      | 4B  | 673,746,810   |
|                      | 4D  | 520,815,567   |
|                      | 5A  | 712,547,961   |
|                      | 5B  | 703,299,309   |
|                      | 5D  | 569,771,178   |
|                      | 6A  | 620,176,429   |
|                      | 6B  | 717,542,863   |
|                      | 7A  | 746,502,734   |
|                      | 7B  | 752,612,656   |
|                      | 7D  | 648,661,963   |
                  








##############################################
# 2 Verify each chr/contig length 
##############################################

mamba activate bio

# 2.1 00_TriticumAestivumRenan
cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/00_TriticumAestivumRenan
FASTAIN="/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/00_TriticumAestivumRenan/Triticum_aestivum_renan.Triticum_aestivum_Renan_v2.1.dna.toplevel.fa.gz"
OUT=/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/00_TriticumAestivumRenan/
seqkit fx2tab -nl ${FASTAIN} -o ${OUT}/Summary_TriticumAestivumRenan.txt
# # chr 21



# 2.2 01_SecaleCereale
cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/01_SecaleCereale
FASTAIN2="/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/01_SecaleCereale/Secale_cereale.Rye_Lo7_2018_v1p1p1.dna.toplevel.fa.gz"
OUT=/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/01_SecaleCereale/
seqkit fx2tab -nl ${FASTAIN2} -o ${OUT}/Summary_Secale_cereale.Rye.txt



# 2.3 02_AvenaSativaOats
cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/02_AvenaSativaOats
FASTAIN3="/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/02_AvenaSativaOats/Avena_sativa_ot3098.Oat_OT3098_v2.dna.toplevel.fa.gz"
OUT=/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/02_AvenaSativaOats/
seqkit fx2tab -nl ${FASTAIN3} -o ${OUT}/Summary_AvenaSativaOats.txt











#############################################################################################################
################  NA ########################################################################################
#############################################################################################################

#----------------------------------------------------------------------
# 1.1 Pinus taeda (Loblolly Pine) - [21 Gb]
#    https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000404065.3/
#----------------------------------------------------------------------
cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/

mkdir 00_PinusTaeda
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/065/GCA_000404065.3_Ptaeda2.0/GCA_000404065.3_Ptaeda2.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/404/065/GCA_000404065.3_Ptaeda2.0/GCA_000404065.3_Ptaeda2.0_genomic.gbff.gz

zcat GCA_000404065.3_Ptaeda2.0_genomic.fna.gz | grep ">" | wc -l  # 1760464

# will need to convert "gbff" into "gtf"
# gffread GCA_000404065.3_Ptaeda2.0_genomic.gbff -T -o Ptaeda2.0.gtf
INFA="/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/z01_PinusTaeda/GCA_000404065.3_Ptaeda2.0_genomic.fna.gz"
OUT=/AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/z01_PinusTaeda/
seqkit fx2tab -nl ${INFA} -o ${OUT}/Summary_Ptaeda2.txt



#------------------------------------------------------------------------------------------
# 1.2 Picea abie [12 Gb]
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/067/695/GCA_900067695.1_Pabies01/
#------------------------------------------------------------------------------------------
cd /AAGI_data/00_nectars_data/A_PROJECTS/sc_ShardSplitRef/05_testing_genomes/

mkdir 01_PiceaAbie
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/067/695/GCA_900067695.1_Pabies01/GCA_900067695.1_Pabies01_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/067/695/GCA_900067695.1_Pabies01/GCA_900067695.1_Pabies01_genomic.gbff.gz 
