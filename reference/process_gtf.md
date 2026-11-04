# Process a GTF file by assigning split-region sequence names

Reads a BED-like "split regions" file and a GTF, assigns each GTF
feature to exactly one split-region interval, and replaces the GTF `Chr`
field with a split-region name of the form `Chr-RegionStart-RegionEnd`.

## Usage

``` r
process_gtf(
  split_regions_bed,
  gtf,
  genome_name,
  genome_version = NULL,
  keep_attributes = NULL,
  out_path = "."
)
```

## Arguments

- split_regions_bed:

  Character: File path to a BED-like file containing split intervals as
  produced by
  [`determine_split_regions()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/determine_split_regions.md).
  The first three columns must be: `Chr`, `RegionStart`, `RegionEnd`.

- gtf:

  Path to a GTF file containing gene annotations. Must contain standard
  9 GTF columns. Comments starting with '#' are ignored.

- genome_name:

  Genome identifier used in the output filename.

- genome_version:

  Genome version string used in the output filename.

- keep_attributes:

  Character vector of attribute keys to keep (e.g.,
  `c("gene_id", "transcript_id", "gene_name")`). If provided, the GTF
  attribute column will be filtered to keep only these fields. Defaults
  to `NULL` (no filtering).

- out_path:

  Directory where the output GTF will be written. Defaults to the
  current working directory (`"."`).

## Value

The function writes the filtered GTF file to local disk at the
user-specified location.

## Details

A feature is considered inside a region when:

- `start >= RegionStart`

- `start < RegionEnd`

- `end <= RegionEnd`

The function aborts if:

- any GTF feature matches more than one region (overlapping regions), or

- any GTF feature matches no region (regions don't fully cover the GTF).

## Output file

The function will automatically name the file based on the input file
and prepends the following string onto the filename:
"B1_FINAL_MODIFIED_GTF\_*genome_name*\_*genome_version*.gtf", placing
the file in the directory specified in the `out_path`.

## Examples

``` r
reg_file <- system.file(
  "extdata/AlgorithmToy_DetermineSplitReg.bed",
  package = "scShardSplitRef",
  mustWork = TRUE
)
gtf_file <- system.file(
  "extdata/AlgorithmToy.gtf",
  package = "scShardSplitRef",
  mustWork = TRUE
)

process_gtf(
  split_regions_bed = reg_file,
  gtf = gtf_file,
  genome_name = "Example",
  genome_version = "v1.0",
  out_path = tempdir()
)
#> ℹ Validating inputs and reading files...
#> All files loaded and validated successfully.
#> ℹ <environment>Finished writing processed GTF: /tmp/RtmpYhLwvz/B1_FINAL_MODIFIED_GTF_Example_v1.0.gtf
```
