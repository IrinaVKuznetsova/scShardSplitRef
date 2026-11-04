# Validate BED and GTF files

Reads and validates BED and GTF files, checking for:

- File existence

- TAB-delimited format with correct number of columns

- Valid BED coordinate ranges (start \< end)

## Usage

``` r
.validate_inputs(bed, gtf)
```

## Arguments

- bed:

  Path to a 3-column BED-like file (chr, start, end). Columns should be
  TAB-delimited with chromosome name, region start (0-based), and region
  end coordinates.

- gtf:

  Path to a GTF file containing gene annotations. Must contain standard
  9 GTF columns. Comments starting with '#' are ignored.

## Value

List with named elements:

- `bed`: Data frame with at least 3 BED columns plus any additional
  columns from the input file.

- `gtf`: Data frame with 9 GTF columns

## Details

Performs four checks in sequence:

1.  Both files exist

2.  GTF has \>= 9 columns

3.  BED has \>= 3 columns

4.  All BED ranges satisfy RegionStart \< RegionEnd

## Examples

``` r
  regions <- system.file("extdata/AlgorithmToy.bed",
                   package = "scShardSplitRef",
                   mustWork = TRUE)

  genes <- system.file("extdata/AlgorithmToy.gtf",
             package = "scShardSplitRef",
             mustWork = TRUE)

  result <- validate_inputs(regions, genes)
#> Error in validate_inputs(regions, genes): could not find function "validate_inputs"

  bed <- result$bed
#> Error: object 'result' not found
  gtf <- result$gtf
#> Error: object 'result' not found

  gtf
#> Error: object 'gtf' not found
  bed
#> Error: object 'bed' not found
```
