# Shift a single boundary past genes

Adjusts one boundary position that falls inside or too close to gene
intervals by shifting it rightward until it satisfies the requested
clearance around genes. Each shift jumps directly to the nearest
shift_by-reachable position that clears the currently-violating
range(s), rather than walking forward one shift_by step at a time, so
the result is identical to a naive per-step walk but computed in
O(number of blocking gene clusters) instead of O(distance / shift_by).

## Usage

``` r
.shift_single_boundary(boundary, gene_ranges, chr_end, shift_by, clearance)
```

## Arguments

- boundary:

  Integer: initial boundary position.

- gene_ranges:

  Data frame with gene ranges (columns 'start', 'end'), already merged
  and sorted ascending by start (and therefore also by end, since ranges
  are non-overlapping).

- chr_end:

  Integer: chromosome end coordinate (upper bound).

- shift_by:

  Integer: step size (bp) for rightward boundary updates.

- clearance:

  Integer: minimum required clearance (bp) from nearby genes on both
  sides.

## Value

Integer: adjusted boundary position. May be the original position if it
did not overlap any gene, or identical to its start position if the
shift would exceed chr_end.

## Details

Candidate gene ranges are also found via two binary searches
([`findInterval()`](https://rdrr.io/r/base/findInterval.html) against
`gene_ranges$start` and `gene_ranges$end`, both sorted ascending since
`gene_ranges` is pre-merged and non-overlapping) rather than a linear
scan of every range on the chromosome, turning an O(n_genes) check into
O(log n_genes + k), where k is the (typically small) number of ranges
actually within `clearance` of `boundary`.

## Examples

``` r
if (FALSE) { # \dontrun{
  adjusted <- .shift_single_boundary(60L, gene_ranges, 300L, 1L, 1L)
} # }
```
