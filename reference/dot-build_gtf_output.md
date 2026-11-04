# Build the output GTF data frame

Selects and renames columns so the output matches GTF layout, replacing
`Chr` with the split-region name stored in `NEW`.

## Usage

``` r
.build_gtf_output(matched)
```

## Arguments

- matched:

  Data frame containing `NEW` and standard GTF columns.

## Value

A data frame with GTF columns where the first column is `Chr`.
