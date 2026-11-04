# scShardSplitRef: Getting Started

## Package Overview

{scShardSplitRef} designed for researchers working with multiome
single-cell data from species that require building a custom reference
genome, particularly when one or more chromosomes/contigs in the
reference FASTA exceed the Cell Ranger ARC limit of
`536.8Mb (2\^29 bases)`. {scShardSplitRef} enables splitting of long
chromosomes/contigs generating shorter fragments that are compatible
with Cell Ranger ARC requirements. The software outputs FASTA and GTF
files containing the split chromosomes/contigs, which can be directly
processed by `cellranger-arc mkref`.

To familiarise yourself with the package and workflow of preparing
reference genome sequence (FASTA file) and gene annotations (GTF file)
for Cell Ranger ARC we suggest starting with the [Guided
Walk-Through](#guided-walk-through). In this tutorial we are using
barley genome (MOREX V3). Example data files are distributed with the
package for example purposes.

## Guided Walk-Through

### Installation

### Prerequisites

| Software        | Version | Purpose                         |
|-----------------|---------|---------------------------------|
| bedtools        | ≥ 2.30  | BED file operations             |
| Cell Ranger ARC | ≥ 2.0   | Single-cell multiome processing |

### Input Files

| File | Format | Description | Comment |
|----|----|----|----|
| Gene annotations | GTF | Gene features and coordinates | Morex v3, subset to 2 chromosomes |
| Centromere coordinates | BED-like | Centromere positions per chromosome (x1) | 2 chromosomes |
| Genome sequence | FASTA | Reference genome sequence | (x2) Not bundled; used for Cell Ranger ARC demo only; download from Ensembl |

> **(x1) Note1:** If centromere coordinates are not available for your
> species, the second column of the BED-like file should be set to `0`.
> If centromere coordinates are available, the second column should be
> set to the centromere coordinate. {scShardSplitRef} will then check
> for any overlaps with genic regions and generate the correct GTF file
> that is compitable with `cellranger-arc mkref`.
>
> **(x2) Note2:** Please note that the FASTA file is not required for
> the package, but is used to demonstrate FASTA preparation for the Cell
> Ranger ARC command in the workflow description. It is not included in
> the repository due to the large file size; download from Ensembl.

  

## Step-by-Step Guide

For this tutorial, we will be using a subset of barley genome (*Hordeum
Vulgare, Morex V3*). We begin by downloading the FASTA and GTF files
from [Ensembl](https://ftp.ensemblgenomes.ebi.ac.uk) and creating a
tab-delimited BED-like file with centromere coordinates from
[GrainGenes](https://graingenes.org/GG3/content/morex-v3-files-2021). We
then assess the formatting of the input GTF and BED files, followed by
checking and, if necessary, modifying the split coordinates in the BED
file. The GTF file is subsequently updated according to the split
position. Finally, the FASTA file is split and formatted for
compatibility with Cell Ranger ARC. All files must be decompressed prior
to use.

  

**What does a GTF file for MOREX V3 look like?**

``` r

# Gene annotations (GTF) file. The format of the GTF file is a tab-delimited 9-column file.
# data.dir = "/inst/extdata/morex3_subset_1H_2H.gtf.gz"

gtf_file <- system.file(
  "extdata/morex3_subset_1H_2H.gtf.gz",
  package = "scShardSplitRef"
)
decompress_gtf <- gzfile(gtf_file, open = "rt")
gtf_df <- read.table(file = decompress_gtf, header = FALSE, sep = "\t")
close(decompress_gtf)
dim(gtf_df)
```

    ## [1] 150638      9

``` r

gtf_df[1:10, ]
```

    ##    V1  V2          V3        V4        V5 V6 V7 V8
    ## 1  2H IPK        gene 601169528 601171444  .  -  .
    ## 2  2H IPK  transcript 601169528 601171444  .  -  .
    ## 3  2H IPK        exon 601169754 601171444  .  -  .
    ## 4  2H IPK         CDS 601169754 601171444  .  -  0
    ## 5  2H IPK start_codon 601171442 601171444  .  -  0
    ## 6  2H IPK        exon 601169528 601169669  .  -  .
    ## 7  2H IPK         CDS 601169531 601169669  .  -  1
    ## 8  2H IPK  stop_codon 601169528 601169530  .  -  0
    ## 9  2H IPK        gene 161203690 161204337  .  -  .
    ## 10 2H IPK  transcript 161203690 161204337  .  -  .
    ##                                                                                                                                                                                                                                                                     V9
    ## 1                                                                                                                                                                                     gene_id HORVU.MOREX.r3.2HG0191020; gene_source IPK; gene_biotype protein_coding;
    ## 2                                                         gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; tag Ensembl_canonical;
    ## 3  gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; exon_number 1; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; exon_id HORVU.MOREX.r3.2HG0191020.1-E1; tag Ensembl_canonical;
    ## 4  gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; exon_number 1; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; protein_id HORVU.MOREX.r3.2HG0191020.1; tag Ensembl_canonical;
    ## 5                                          gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; exon_number 1; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; tag Ensembl_canonical;
    ## 6  gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; exon_number 2; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; exon_id HORVU.MOREX.r3.2HG0191020.1-E2; tag Ensembl_canonical;
    ## 7  gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; exon_number 2; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; protein_id HORVU.MOREX.r3.2HG0191020.1; tag Ensembl_canonical;
    ## 8                                          gene_id HORVU.MOREX.r3.2HG0191020; transcript_id HORVU.MOREX.r3.2HG0191020.1; exon_number 2; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; tag Ensembl_canonical;
    ## 9                                                                                                                                                                                     gene_id HORVU.MOREX.r3.2HG0134140; gene_source IPK; gene_biotype protein_coding;
    ## 10                                                        gene_id HORVU.MOREX.r3.2HG0134140; transcript_id HORVU.MOREX.r3.2HG0134140.1; gene_source IPK; gene_biotype protein_coding; transcript_source IPK; transcript_biotype protein_coding; tag Ensembl_canonical;

**What does a BED-like centromere file for MOREX V3 look like?**

``` r

# The format of the centromere coordinates file is a tab-delimited 3-column file, where
#   the 1st column contains chromosome/contig names,
#   the 2nd column contains the centromere coordinate (or 0 if this information is not available for your species),
#   the 3rd column contains the chromosome/contig end coordinates.
# Manually formatted centromere file (please check the expected format below):
# data.dir = "/inst/extdata/morex3_subset_1H_2H.bed"

bed_file <- system.file(
  "extdata/morex3_subset_1H_2H.bed",
  package = "scShardSplitRef"
)
centr_bed_df <- read.table(file = bed_file, header = FALSE, sep = "\t")
dim(centr_bed_df)
```

    ## [1] 2 3

``` r

centr_bed_df
```

    ##   V1        V2        V3
    ## 1 1H 206486643 516505932
    ## 2 2H 301293086 665585731

### Step-1

The goal of the
[`determine_split_regions()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/determine_split_regions.md)
function is to create appropriate split regions. It initially checks GTF
and BED-like input files, then prepares split coordinates in the
tab-delimited BED-like format. This function ensures that chromosomes
are not split in the middle of a gene region. The BED-like file with
centromere coordinates if available, or zeros if not. This file is
manually created. All files need to be decompressed.

### The Help Page for the `determine_split_regions()` Function

``` r

?determine_split_regions
```

### Check formatting of input files and build BED-like split intervals

``` r

# Load the library
library(scShardSplitRef)

# Load MOREX V3 (subset to 1H, 2H)
# Note1: all files need to be decompressed for use in scShardSplitRef and Cell Ranger ARC
# Note2: use correctly formatted BED file with centromere coordinates
# Task: Check formatting of input files. Build BED-like split intervals using genomic regions (BED) and gene annotations (GTF) input files

# Read GTF file
gtf_path <- system.file(
  "extdata/morex3_subset_1H_2H.gtf.gz",
  package = "scShardSplitRef"
)
tmp_gtf_path <- tempfile(fileext = ".gtf")
con_in <- gzcon(file(gtf_path, open = "rb"))
con_out <- file(tmp_gtf_path, open = "wb")
writeBin(readBin(con_in, raw(), n = 1e8), con_out)
close(con_in)
close(con_out)

# Read BED file
bed_path <- system.file(
  "extdata/morex3_subset_1H_2H.bed",
  package = "scShardSplitRef"
)

# Define output path
output_file <- file.path(tempdir(), "split")
dir.create(output_file, recursive = TRUE, showWarnings = FALSE)
output_file_name_dir <- file.path(
  output_file,
  "morex3_subset_1H_2H_Split_Regions.bed"
)

# Determine split regions
get_split_reg <- determine_split_regions(
  bed = bed_path,
  gtf = tmp_gtf_path,
  output_bed = output_file_name_dir,
  limit = 2L^29L,
  shift_by = 20L,
  clearance = 20L
)
```

    ## Preparing split regions from
    ## /home/runner/work/_temp/Library/scShardSplitRef/extdata/morex3_subset_1H_2H.bed
    ## and /tmp/RtmppT2HFo/file1bbc14d2f919.gtf into
    ## /tmp/RtmppT2HFo/split/morex3_subset_1H_2H_Split_Regions.bed (`limit` =
    ## 536870912, `shift_by` = 20, `clearance` = 20).
    ## Built 4 split intervals across 2 chromosomes and wrote them to
    ## /tmp/RtmppT2HFo/split/morex3_subset_1H_2H_Split_Regions.bed.

``` r

get_split_reg
```

    ## [1] "/tmp/RtmppT2HFo/split/morex3_subset_1H_2H_Split_Regions.bed"

**What does the generated BED file look like?\`**

``` r

bed_f <- system.file(
  "extdata/morex3_subset_1H_2H_Split_Regions.bed",
  package = "scShardSplitRef"
)
split_region_bed_df <- read.table(file = bed_f, header = FALSE, sep = "\t")
dim(split_region_bed_df)
```

    ## [1] 4 3

``` r

split_region_bed_df
```

    ##   V1        V2        V3
    ## 1 1H         0 206486643
    ## 2 1H 206486643 516505932
    ## 3 2H         0 301293086
    ## 4 2H 301293086 665585731

### Step-2

Generate a gene annotation (GTF) file in compatible to Cell Ranger ARC
format. The first column of the GTF file should be in the following form
`Chr-RegionStart-RegionEnd`. Note that the output file name is hardcoded
and starts with
“B1_FINAL_MODIFIED_GTF\_\[*genome_name*\]\_\[*genome_version*\].gtf”.

#### The Help Page for the `process_gtf()` Function

``` r

?process_gtf
```

#### Generate Gene Annotation (GTF) File in Compatible to Cell Ranger ARC Format

``` r

# Process GTF file and write file to tempdir()
print(output_file_name_dir)
```

    ## [1] "/tmp/RtmppT2HFo/split/morex3_subset_1H_2H_Split_Regions.bed"

``` r

print(tmp_gtf_path)
```

    ## [1] "/tmp/RtmppT2HFo/file1bbc14d2f919.gtf"

``` r

process_gtf(
  split_regions_bed = output_file_name_dir,
  gtf = tmp_gtf_path,
  genome_name = "subset_barley",
  genome_version = "3",
  out_path = tempdir()
)
```

    ## ℹ Validating inputs and reading files...

    ## All files loaded and validated successfully.
    ## ℹ <environment>Finished writing processed GTF: /tmp/RtmppT2HFo/B1_FINAL_MODIFIED_GTF_subset_barley_3.gtf

> The output file is in the GTF format. Note that the first column has
> been reformatted, the regions are now split.  
> The chromosome name is formatted as `1H-0-206486643`, which is the
> region from the beginning of the chromosome to the split position;  
> `1H-206486643-516505932` is the region from the split position to the
> end of the chromosome.

**What does the generated GTF file look like?\`**

``` r

gtf_modified <- list.files(
  tempdir(),
  pattern = "B1_FINAL_MODIFIED_GTF_subset.*\\.gtf$",
  full.names = TRUE
)

processed_gtf <- read.table(file = gtf_modified, header = FALSE, sep = "\t")
print(gtf_modified)
```

    ## [1] "/tmp/RtmppT2HFo/B1_FINAL_MODIFIED_GTF_subset_barley_3.gtf"

``` r

dim(processed_gtf)
```

    ## [1] 150638      9

``` r

sep <- as.data.frame(matrix("...", nrow = 1, ncol = 5))
colnames(sep) <- colnames(processed_gtf)[1:5]

kable(rbind(
  processed_gtf[200:203, 1:5],
  sep,
  processed_gtf[(nrow(processed_gtf) - 2):nrow(processed_gtf), 1:5]
))
```

|        | V1                     | V2  | V3              | V4        | V5        |
|:-------|:-----------------------|:----|:----------------|:----------|:----------|
| 200    | 2H-0-301293086         | IPK | CDS             | 146924435 | 146924961 |
| 201    | 2H-0-301293086         | IPK | stop_codon      | 146924962 | 146924964 |
| 202    | 2H-301293086-665585731 | IPK | gene            | 130299937 | 130300482 |
| 203    | 2H-301293086-665585731 | IPK | transcript      | 130299937 | 130300482 |
| 1      | …                      | …   | …               | …         | …         |
| 150636 | 1H-206486643-516505932 | IPK | stop_codon      | 291315417 | 291315419 |
| 150637 | 1H-206486643-516505932 | IPK | five_prime_utr  | 291310927 | 291310985 |
| 150638 | 1H-206486643-516505932 | IPK | three_prime_utr | 291315420 | 291315710 |

  
**{scShardSplitRef} generated a GTF file with split coordinates that can
be used in Cell Ranger ARC.**  

### Step-3

Here we show the workflow on obtaining FASTA files from the Ensembl FTP
server, explaining how to split FASTA and keep the header naming
convention in the same style as the GTF file generated by the
`scShardSplitRef` package.

> **Important:** `bedtools getfasta` writes FASTA headers as
> `Chr:RegionStart-RegionEnd` by default. You must convert `:` to `-` so
> FASTA sequence names match the GTF sequence names generated by
> [scShardSplitRef](https://irinavkuznetsova.github.io/scShardSplitRef/).
> If this is skipped, `cellranger-arc mkref` may fail with sequence-name
> mismatch errors.

``` bash
# Get FASTA from Ensembl
# Hordeum vulgare, (release-62), Ensembl
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.1H.fa.gz
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/hordeum_vulgare/dna/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.2H.fa.gz
```

``` bash
# Decompress FASTA
gunzip Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.1H.fa.gz
gunzip Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.2H.fa.gz
```

``` bash
# Concatenate into one file
cat Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.1H.fa Hordeum_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.2H.fa >morex3_subset_1H_2H.fasta
```

``` bash
# Please note that due to the file size limitations, the barley genome FASTA file is not available on our GitHub repository at /inst/extdata/
IN_FA="/data/morex3_subset_1H_2H.fasta"
BED_FIN_OUT="/inst/extdata/morex3_subset_1H_2H_Split_Regions.bed"
FASTA_FIN_OUT="/data/"

# Split FASTA
bedtools getfasta -fi ${IN_FA} -bed ${BED_FIN_OUT} -fo ${FASTA_FIN_OUT}/morex3_SPLIT_1H_2H.fasta

# CRITICAL: Do not skip this rename step; Cell Ranger ARC expects FASTA names to match GTF names exactly
# Convert FASTA headers from `Chr:RegionStart-RegionEnd` to `Chr-RegionStart-RegionEnd`
sed '/^>/s/:/-/g' ${FASTA_FIN_OUT}/morex3_SPLIT_1H_2H.fasta >${FASTA_FIN_OUT}/morex3_split_FINAL_1H_2H.fasta

# Validate
seqkit fx2tab --length --name ${FASTA_FIN_OUT}/morex3_split_FINAL_1H_2H.fasta

1H-0-206486643 206486643
1H-206486643-516505932 310019289
2H-0-301293086 301293086
2H-301293086-665585731 364292645
```

**What does the generated FASTA file look like?\`**

    # morex3_split_FINAL_1H_2H.fasta

    >1H-0-206486643
    GCATTTGGTTGCTTTTATTCTAAACTA...
    >1H-206486643-516505932
    ATCAGTACATGAGCTAGTATTAGAGGCTAGACTAGGGATCAATTGGTGTTTGTA...
    >2H-0-301293086
    GCCCTGAC...

### Step-4

The final step demonstrates all required input files for the Cell Ranger
ARC command. The genome sequence (FASTA) and gene annotation (GTF) files
are now split and formatted for Cell Ranger ARC, with file paths
specified in the .config file, and Cell Ranger ARC can now be run.

#### Prepare `config` File

``` bash
cat 03_cellranger_arc_barley.config
```

``` json
{
    "organism": "barley",
    "genome": ["MOREX3_subset_v3"],
    "input_fasta": ["/data/morex3_split_FINAL_1H_2H.fasta"],
    "input_gtf": ["/extdata/B1_FINAL_MODIFIED_GTF_subset_1H_2H_morex_3.gtf"]
}
```

> **Note:** In the full genome (not the subset), you need to add
> information about mitochondrial and chloroplast to the
> “non_nuclear_contigs” field in the config file:
> `"non_nuclear_contigs": ["Mt-0-525599", "Pt-0-115974"]` in the config
> file

#### Run Cell Ranger ARC

``` bash
CONFIG_IN="03_cellranger_arc_barley.config"

cellranger-arc mkref \
  --config=${CONFIG_IN} \
  --nthreads=24 \
  >mkref.log 2>&1
```

## Citation

Please cite this tool if you use it

``` r

citation("scShardSplitRef")
```

    ## To cite package 'scShardSplitRef' in publications use:
    ## 
    ##   Kuznetsova I, Pembleton L, Curci L, Sparks A (2026).
    ##   _scShardSplitRef: Reference Genome Splitting for Single-Cell Multiome
    ##   Data_. R package version 0.1.0,
    ##   <https://irinavkuznetsova.github.io/scShardSplitRef/>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {scShardSplitRef: Reference Genome Splitting for Single-Cell Multiome Data},
    ##     author = {Irina Kuznetsova and Luke Pembleton and Luca Curci and Adam H. Sparks},
    ##     year = {2026},
    ##     note = {R package version 0.1.0},
    ##     url = {https://irinavkuznetsova.github.io/scShardSplitRef/},
    ##   }

## License

This project is licensed under the GPL-3 license

## Contributing

Please refer to the CONTRIBUTING.md If you’d like to help improve this
tool, please send your suggestions via email. Please clearly outline the
feature(s) you think should be added. As the tool develops, new features
can be added as separate package functions and contributed back to the
project.  
Remember that improving tools is important, but keeping them neat and
not over-complicated is just as essential.

## Acknowledgements

GRDC funded this project through the Analytics for the Australian Grains
Industry (AAGI) investment, CUR2210-005OPX (AAGI-CU).

This work was supported by resources provided by the [Pawsey
Supercomputing Research Centres](https://pawsey.org.au/) [Nimbus
Research Cloud](https://doi.org/10.48569/v0j3-qd51), with funding from
the Australian Government and the Government of Western Australia.

This research used the Australian Research Data Commons (ARDC) ARDC
Nectar Research Cloud Service [ARDC Nectar Research Cloud \| Australian
Research Data
Commons](https://ardc.edu.au/services/ardc-nectar-research-cloud/). The
ARDC is enabled by the Australian Government’s National Collaborative
Research Infrastructure Strategy (NCRIS).
