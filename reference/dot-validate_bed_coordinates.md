# Validate BED coordinate format

Checks that start and end coordinates satisfy BED format requirements:
start \>= 0 and start \< end.

## Usage

``` r
.validate_bed_coordinates(start, end, chr)
```

## Arguments

- start:

  Integer: region start position (0-based).

- end:

  Integer: region end position (exclusive).

- chr:

  Character: chromosome name for error messages.

## Value

Invisible TRUE if valid. Raises error if invalid.

## Examples

``` r
if (FALSE) { # \dontrun{
  .validate_bed_coordinates(0L, 100L, "chr1")
} # }
```
