# Identify rows where a feature lies fully within a split region

Applies the matching rule used by
[`process_gtf()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/process_gtf.md).
A feature `start`, `end` should be fully inside a region: `RegionStart`,
`RegionEnd`: `start >= RegionStart` `end <= RegionEnd`

## Usage

``` r
.rows_within_region(df)
```

## Arguments

- df:

  Data frame containing `start`, `end`, `RegionStart`, `RegionEnd`.

## Value

Logical vector of length `nrow(df)`.

## Details

Retained for reference/tests; the main pipeline now enforces this rule
directly inside
[`.match_features_to_regions()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-match_features_to_regions.md)
without a full cross join.
