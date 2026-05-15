---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->




## **scShardSplitRef**  
*scShardSplitRef* designed for researchers working with multiome single-cell data from species that require building a custom reference genome, particularly when one or more chromosomes/contigs in the reference FASTA exceed the Cell Ranger ARC limit of 536.8Mb (2\^29 bases). *scShardSplitRef* enables splitting of long chromosomes/contigs generating shorter fragments that are compitable with Cell Ranger ARC requirements. The software outputs FASTA and GTF files containing the split chromosomes/contigs, which can be directly processed by `cellranger-arc mkref`.

To familiarise yourself with the package and workflow of preparing reference genome sequence (FASTA) and gene annotations (GTF) files for Cell Ranger ARC we suggest starting with a guided walk through. 


## **Problem overview**  
*Cell Ranger* provides pre-built reference genomes for common species such as human and mouse. For less common species, users must build a custom reference. For example, for multiome ATAC + Gene Expression sequencing data the command is `cellranger-arc mkref`tha requires geneome sequence (FASTA) and gene annotations (GTF) files. When working with species that have very large chromosomes or scaffolds, users may encounter the following error:  
>*Due to limitations of the BAM index format, a contig in the reference FASTA file cannot exceed *536.8Mb (2^29 bases)*. If a contig exceeds that size you will have to split it into smaller contigs and make corresponding modifications to the GTF file.* 
[Source: 10x Genomics Documentation](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/inputs/mkref)  

>![screenshot of the error](man/figures/scr_10x_CR_arc_error.png)  


Barley genome has several chromosomes (*2H=665.59Mb, 3H=621.52Mb, 4H=610.33Mb, 5H=588.22Mb, 6H=561.79Mb, 7H=632.54Mb*) that exceed this size limit. As a result, `cellranger-arc mkref` cannot build a reference unless "oversized" chromosomes are divided into smaller fragments.  




## scShardSplitRef guided walk through

#### Installation  
```{r}
# Install the package
# install.packages("pak")

# Load the package 
library("dplyr")
library("scShardSplitRef")
```

For this tutorial, we will be using a genome sequence (FASTA file) and gene annotations (GTF file) of the barley genome (*Hordeum Vulgare, Morex V3*) available from [Ensembl](https://ftp.ensemblgenomes.ebi.ac.uk) and centromere coordinates from [GrainGenes](https://graingenes.org/GG3/content/morex-v3-files-2021).

We begin by assessing the formatting of the input GTF and BED files provided by the user, followed by checking and, if necessary, modifying the split coordinates in the BED file. The GTF file is then updated according to the split position. Finally, FASTA file is split and formatted for compatibility with Cell Ranger ARC.
<details>
<summary><b>What do GTF and centromere files for MOREX V3 look like?</b></summary>

```{r, eval = FALSE}
# Gene annotations (GTF) file: 
# data.dir = "/inst/extdata/01_Hordeum_vulgare.Mt_Pt_v49_Barley_v62_COMBINED.gtf"

gtf_df <- read.table(file = "01_Hordeum_vulgare.Mt_Pt_v49_Barley_v62_COMBINED.gtf",
  header = FALSE,
  sep = "\t")
dim(gtf_df)     ### 524562 x 9
gtf_df[1:10, ]
```
The GTF file conatins 9 columns.
```
#!genome-build MorexV3_pseudomolecules_assembly
#!genome-version MorexV3_pseudomolecules_assembly
#!genome-date 2021-04
#!genome-build-accession GCA_904849725.1
#!genebuild-last-updated 2021-04       
2H      IPK     gene    601169528       601171444       .       -       .       gene_id "HORVU.MOREX.r3.2HG0191020"; gene_source "IPK"; gene_biotype "protein_coding";
2H      IPK     transcript      601169528       601171444       .       -       .       gene_id "HORVU.MOREX.r3.2HG0191020"; transcript_id "HORVU.MOREX.r3.2HG0191020.1"; gene_source "IPK"; gene_biotype "protein_coding"; transcript_source "IPK"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
2H      IPK     exon    601169754       601171444       .       -       .       gene_id "HORVU.MOREX.r3.2HG0191020"; transcript_id "HORVU.MOREX.r3.2HG0191020.1"; exon_number "1"; gene_source "IPK"; gene_biotype "protein_coding"; transcript_source "IPK"; transcript_biotype "protein_coding"; exon_id "HORVU.MOREX.r3.2HG0191020.1-E1"; tag "Ensembl_canonical";
2H      IPK     CDS     601169754       601171444       .       -       0       gene_id "HORVU.MOREX.r3.2HG0191020"; transcript_id "HORVU.MOREX.r3.2HG0191020.1"; exon_number "1"; gene_source "IPK"; gene_biotype "protein_coding"; transcript_source "IPK"; transcript_biotype "protein_coding"; protein_id "HORVU.MOREX.r3.2HG0191020.1"; tag "Ensembl_canonical";
2H      IPK     start_codon     601171442       601171444       .       -       0       gene_id "HORVU.MOREX.r3.2HG0191020"; transcript_id "HORVU.MOREX.r3.2HG0191020.1"; exon_number "1"; gene_source "IPK"; gene_biotype "protein_coding"; transcript_source "IPK"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
```



```{r}
# Manually formatted centromere file (please check what it should look like below and check **Step-0** description): 
# data.dir = "/inst/extdata/02_MorexV3_centromere_coordinates.bed"

centr_bed_df <- read.table(file = "02_MorexV3_centromere_coordinates.bed",
  header = FALSE,
  sep = "\t")
dim(centr_bed_df)  ## 292 x 3
centr_bed_df[1:10, ]
```

```         
1H      206486643       516505932
2H      301293086       665585731
3H      267852507       621516506
4H      276149121       610333535
5H      204878572       588218686
6H      256319444       561794515
7H      328847192       632540561
CAJHDD010000004.1       0       972986
CAJHDD010000001.1       0       680018
CAJHDD010000005.1       0       517010
```

</details>



#### Step-0
All files need to be decompressed. The format of the centromere coordinate file is a tab-delimited 3-column file, where the 1st column contains chromosome/contig names, the 2nd column contains the centromere coordinate, and the 3rd column contains chromosome/contig region end coordinates.

```{r}
# Load MOREX V3 files
# Note1: all files need to be decompressed
# Note2: use correctly formatted BED file with centromere coordinates 
GTF_IN = "/inst/extdata/01_Hordeum_vulgare.Mt_Pt_v49_Barley_v62_COMBINED.gtf"
BED_CENTROM_IN = "/inst/extdata/02_MorexV3_centromere_coordinates.bed"
```

#### Step-1
The initial step is to check if required files are provided in the correct format.

```{r}
?validate_inputs
input_verif <- validate_inputs(bed = BED_CENTROM_IN, 
                               gtf = GTF_IN)
```
If GTF and BED formats are correct, a confirmation message will be printed to the console.
```         
All files loaded and validated successfully.
```

#### Step-2
The next step is to prepare split coordinates in BED-like format. This function ensures that chromosomes are not split in the middle of a gene region. To achieve this, we first check whether the centromere coordinates overlap with any gene features, and if they do, we shift the split position by a specified length in base pairs (bp).

```{r}
# Build BED-like split intervals from genomic regions (BED) and gene annotations (GTF)
OUTPUT_DIR_BED = "/outdata"
?determine_split_regions
get_split_reg <- determine_split_regions(bed=BED_CENTROM_IN, 
                        gtf = GTF_IN, 
                        file.path("/outdata/", "B1_Split_regions.bed"), 
                        limit = 2L^29L, 
                        shift_by = 20L)
```

<details>
<summary><b>What do the results "B1_Split_regions.bed" looks like?</b></summary>

```         
1H  0   206486643
1H  206486643   516505932
2H  0   301293086
2H  301293086   665585731
3H  0   267852507
3H  267852507   621516506
4H  0   276149121
4H  276149121   610333535
5H  0   204878572
5H  204878572   588218686
6H  0   256319444
6H  256319444   561794515
7H  0   328847192
7H  328847192   632540561
CAJHDD010000004.1   0   972986
CAJHDD010000001.1   0   680018
```

</details>

#### Step-3
Generate gene annotation (GTF) file in compatible to Cell Ranger ARC format. The first column of the GTF file should be in the following form `Chr-RegionStart-RegionEnd`. Note that the output file name is hardcoded and starts with "B1_FINAL_MODIFIED_GTF_[*genome_name*]_[*genome_version*].gtf". 

```{r}
GTF_FIN_OUT = "/outdata/"
GTF_processed <- process_gtf(split_regions_bed = "/outd ata/B1_Split_regions.bed", 
                          gtf = GTF_IN, 
                          genome_name = "barley", 
                          genome_version = "3", 
                          out_path = GTF_FIN_OUT )
```

<details>
<summary><b>What do the results in "B1_FINAL_MODIFIED_GTF_barley_3.gtf" looks like?</b></summary>
The output file is in GTF format. Note that the first column has been reformatted — the regions are now split. The chromosome name is formatted as ```1H-0-206486643```, which is the region from the beginning of the chromosome to the split position; ```1H-206486643-516505932`` is the region from the split position to the end of the chromosome. 
  
```         
1H-0-206486643	IPK	gene	151788036	151791331	.	-	.	gene_id "HORVU.MOREX.r3.1HG0030850"; gene_source "IPK"; gene_biotype "protein_coding";
1H-0-206486643	IPK	transcript	151788036	151791331	.	-	.	gene_id "HORVU.MOREX.r3.1HG0030850"; transcript_id "HORVU.MOREX.r3.1HG0030850.1"; gene_source "IPK"; gene_biotype "protein_coding"; transcript_source "IPK"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
...
1H-206486643-516505932	IPK	gene	132967560	132971272	.	-	.	gene_id "HORVU.MOREX.r3.1HG0050880"; gene_source "IPK"; gene_biotype "protein_coding";
1H-206486643-516505932	IPK	transcript	132967560	132971272	.	-	.	gene_id "HORVU.MOREX.r3.1HG0050880"; transcript_id "HORVU.MOREX.r3.1HG0050880.1"; gene_source "IPK"; gene_biotype "protein_coding"; transcript_source "IPK"; transcript_biotype "protein_coding"; tag "Ensembl_canonical";
...
```
</details>


#### Step-4
Final steps is to split FASTA file, the fatsa headers should have the same chrosomome annotations as GTF file. 

```{linux}
# Please note that due to file size limitations, the barley genome FASTA file is not available on our GitHub repository at /inst/extdata/
IN_FA = "01_Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly_COMBINED_Mt_Pt_v49_Barley_v62.fa"
BED_FIN_OUT = "/outd ata/B1_Split_regions.bed"
FASTA_FIN_OUT = "/outdata/"

# Split FASTA 
bedtools getfasta -fi ${IN_FA} -bed ${BED_FIN_OUT} -fo ${FASTA_FIN_OUT}/B1_FASTA_split.fasta

# Convert FASTA heading to Chr-RegionStart-RegionEnd form
sed '/^>/s/:/-/g' ${FASTA_FIN_OUT}/B1_FASTA_split.fasta > ${FASTA_FIN_OUT}/B1_FASTA_split_final.fasta

```

<details>

<summary><b>What do the results in "B1_FASTA_split_final.fasta" looks like?</b></summary>

```         
>1H-0-206486643
ATTTGTAGTGCTTTTCAATTTCAGGGTCAA.....
>1H-206486643-516505932
ACTGGCCAAAATAGATCAAAATTGCGAGTTTTGACGAGTTCCCCGTAAGCGGACTTCGG...........
>2H-0-301293086
GTGAGAATAGGCC.........
>2H-301293086-665585731
...
```

</details>

#### Step-5
The genome sequence (FASTA) and gene annotation (GTF) files are now split and formatted for Cell Ranger ARC. The file paths are specified in the ```.config`` file and Cell Ranger ARC can now be run.

<details>
<summary><b>What do the ".config" file looks like?</b></summary>

  ```
      
{
    organism: "barley"
    genome: ["MOREX3_v3"]
    input_fasta: ["/outdata/B1_FASTA_split_final.fasta"]
    input_gtf: ["/outdata/B1_FINAL_MODIFIED_GTF_barley_3.gtf"]
    non_nuclear_contigs: ["Mt-0-525599", "Pt-0-115974"]
}

```
 </details> 



## R package session information 

<details>
<summary><b>sessionInfo()</b></summary>
  
```
R version 4.5.3 (2026-03-11)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 22.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_AU.UTF-8       LC_NUMERIC=C               LC_TIME=en_AU.UTF-8        LC_COLLATE=en_AU.UTF-8    
 [5] LC_MONETARY=en_AU.UTF-8    LC_MESSAGES=en_AU.UTF-8    LC_PAPER=en_AU.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_AU.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] scShardSplitRef_0.0.0.9000

loaded via a namespace (and not attached):
 [1] xml2_1.5.2        stringi_1.8.7     digest_0.6.39     magrittr_2.0.5    evaluate_1.0.5    pkgload_1.5.2    
 [7] fastmap_1.2.0     rprojroot_2.1.1   processx_3.9.0    pkgbuild_1.4.8    sessioninfo_1.2.3 brio_1.1.5       
[13] urlchecker_1.0.1  ps_1.9.3          promises_1.5.0    purrr_1.2.2       codetools_0.2-18  cli_3.6.6        
[19] shiny_1.13.0      rlang_1.2.0       pak_0.9.5         commonmark_2.0.0  ellipsis_0.3.3    remotes_2.5.0    
[25] withr_3.0.2       cachem_1.1.0      yaml_2.3.12       devtools_2.4.5    otel_0.2.0        devtag_0.0.0.9000
[31] tools_4.5.3       roxyglobals_1.0.0 memoise_2.0.1     httpuv_1.6.17     vctrs_0.7.3       R6_2.6.1         
[37] mime_0.13         lifecycle_1.0.5   fs_1.6.3          htmlwidgets_1.6.4 usethis_3.2.1     miniUI_0.1.2     
[43] pkgconfig_2.0.3   desc_1.4.3        callr_3.7.6       pillar_1.11.1     later_1.4.8       glue_1.8.1       
[49] profvis_0.4.0     Rcpp_1.1.1-1.1    xfun_0.57         tibble_3.3.1      rstudioapi_0.18.0 knitr_1.51       
[55] xtable_1.8-8      htmltools_0.5.9   rmarkdown_2.31    testthat_3.3.2    compiler_4.5.3    roxygen2_8.0.0 

```

 </details> 



## Citation
Please cite this tool if you use it

``` {r}
library("scShardSplitRef")
citation("scShardSplitRef")
#> Warning in citation("scShardSplitRef"): could not determine year for 'scShardSplitRef'
#> from package DESCRIPTION file
#> To cite package 'scShardSplitRef' in publications use:
#> 
#>   Kuznetsova I, Pembleton L, Curci L, Sparks A (????). _scShardSplitRef: Build
#>   Custom Reference Genome for Single-Cell Multiome ATAC and Gene Expression Data
#>   of Large-genome Species._. R package version 0.0.0.9000.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {scShardSplitRef: Build Custom Reference Genome for Single-Cell Multiome ATAC and
#> Gene Expression Data of Large-genome Species.},
#>     author = {Irina Kuznetsova and Luke Pembleton and Luca Curci and Adam H. Sparks},
#>     note = {R package version 0.0.0.9000},
#>   }
```


## License
This project is licensed under the [GPL-3.0 license](LICENSE.md).  



## Contributing  
If you'd like to help improve this tool, please send your suggestions via email.  
Please clearly outline the feature(s) you think should be added. As the tool develops, new features can be added as separate package functions and contributed back to the project.    
Remember that improving tools is important, but keeping them neat and not over-complicated is just as essential. 


## Acknowledgements  
This work was supported by resources provided by the [Pawsey Supercomputing Research Centres](https://pawsey.org.au/) Nimbus Research Cloud (https://doi.org/10.48569/v0j3-qd51), with funding from the Australian Government and the Government of Western Australia.

GRDC (?)

