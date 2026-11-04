# Validate BED coordinate ranges

Checks that all BED regions satisfy the standard requirement that
RegionStart \< RegionEnd (0-based half-open intervals).

## Usage

``` r
.validate_bed_ranges(bed)
```

## Arguments

- bed:

  Data frame with RegionStart and RegionEnd columns.

## Value

Invisible TRUE if all ranges are valid. Raises error with problematic
row numbers otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
  .validate_bed_ranges(bed_df)
} # }
```
