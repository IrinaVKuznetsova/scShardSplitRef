# Prepare split regions by adding `NEW` region names

Adds a `NEW` column to the split-regions data frame, formatted as
`Chr-RegionStart-RegionEnd`.

## Usage

``` r
.prepare_split_regions(regions)
```

## Arguments

- regions:

  Data frame with columns `Chr`, `RegionStart`, `RegionEnd`.

## Value

The input data frame with an added `NEW` column.
