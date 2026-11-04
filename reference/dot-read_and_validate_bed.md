# Read, validate, and parse a BED file in a single pass

Read, validate, and parse a BED file in a single pass

## Usage

``` r
.read_and_validate_bed(path)
```

## Arguments

- path:

  Character: path to BED file.

## Value

Data frame with first 3 columns named Chr, RegionStart, RegionEnd
(integer), plus any additional columns from the input file.

## Examples

``` r
if (FALSE) { # \dontrun{
  bed_df <- .read_and_validate_bed("regions.bed")
} # }
```
