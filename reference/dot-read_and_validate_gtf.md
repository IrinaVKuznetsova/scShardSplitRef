# Read, validate, and parse a GTF file in a single pass

Read, validate, and parse a GTF file in a single pass

## Usage

``` r
.read_and_validate_gtf(path)
```

## Arguments

- path:

  Character: path to GTF file.

## Value

Data frame with 9 columns named: Chr, source, feature, start, end,
score, strand, frame, attribute (start/end coerced to integer).

## Examples

``` r
if (FALSE) { # \dontrun{
  gtf_df <- .read_and_validate_gtf("annotations.gtf")
} # }
```
