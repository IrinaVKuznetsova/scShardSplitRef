# Assert that the size limit is feasible

Checks that the requested size `limit` is larger than any single gene
with two-sided boundary clearance. Raises an informative error if not
feasible.

## Usage

``` r
.assert_feasible_limit(genes, limit, clearance)
```

## Arguments

- genes:

  Data frame with gene annotations (columns 'chr', 'start', 'end').

- limit:

  Integer: requested maximum interval width.

- clearance:

  Integer: minimum required clearance (bp) from nearby genes on both
  sides.

## Value

Invisible TRUE if feasible. Raises error with remediation advice if not.

## Details

Merges overlapping genes per chromosome via
[`.merge_all_gene_ranges()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-merge_all_gene_ranges.md),
then checks if any merged range plus two-sided clearance (2 \*
`clearance`) exceeds the limit. Reports the worst (largest) blocking
interval and suggests an appropriate minimum limit value.

## Examples

``` r
if (FALSE) { # \dontrun{
  .assert_feasible_limit(genes, limit = 50000, clearance = 1)
} # }
```
