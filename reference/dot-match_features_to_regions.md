# Match GTF features to split regions

Assigns each GTF feature to the single split-region interval that fully
contains it, per chromosome.

## Usage

``` r
.match_features_to_regions(gtf, regions)
```

## Arguments

- gtf:

  Data frame with columns `Chr`, `start`, `end`, `unique_ID` (plus other
  standard GTF columns).

- regions:

  Data frame with columns `Chr`, `RegionStart`, `RegionEnd`, `NEW`.

## Value

`gtf` with three additional columns, `NEW`, `RegionStart`, `RegionEnd`,
populated for matched rows and `NA` for unmatched rows. Row order is
identical to the input `gtf`.

## Details

Rather than joining every feature against every split region on its
chromosome (a full cross join, which duplicates each feature once per
shard on that chromosome before being filtered back down – expensive for
genome-scale GTFs split into many shards), this does a direct,
vectorized per-chromosome lookup: because split regions are asserted to
be non-overlapping (see
[`.regions_overlap()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-regions_overlap.md)),
each feature start position can belong to at most one region, found via
a single sorted-boundary
[`findInterval()`](https://rdrr.io/r/base/findInterval.html) call
instead of a full join. This turns an O(features \*
regions_per_chromosome) join into O(features log
regions_per_chromosome).

A feature is matched to a region only when it is fully contained
(`start >= RegionStart` and `end <= RegionEnd`); a feature whose start
falls in one region but whose end crosses into the next (or that falls
entirely before/after all regions on its chromosome) is left unmatched,
exactly as under the original cross-join-and-filter approach.
