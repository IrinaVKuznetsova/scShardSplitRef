# Check whether any split regions overlap within a chromosome

The fast interval-assignment used to match GTF features to split regions
(see
[`.match_features_to_regions()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-match_features_to_regions.md))
assumes at most one split region can contain any given position on a
chromosome. This checks that invariant upfront, cheaply (region count is
typically small relative to feature count), so
[`process_gtf()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/process_gtf.md)
can choose the fast path when it holds and fall back to the original,
slower-but-exhaustive matching when it doesn't – which is also the only
case where the ambiguous-match error actually fires, so falling back
there costs nothing in practice.

## Usage

``` r
.regions_overlap(regions)
```

## Arguments

- regions:

  Data frame with columns `Chr`, `RegionStart`, `RegionEnd`.

## Value

Logical scalar: `TRUE` if any chromosome has overlapping region
intervals, `FALSE` otherwise.
