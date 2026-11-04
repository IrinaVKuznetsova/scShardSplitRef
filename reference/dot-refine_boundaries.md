# Iteratively refine boundaries

Core algorithm that repeatedly adjusts boundaries by shifting past genes
and subdividing oversized intervals until convergence.

## Usage

``` r
.refine_boundaries(
  boundaries,
  gene_ranges,
  chr_end,
  limit,
  shift_by,
  clearance = 1L,
  max_iter = 100L
)
```

## Arguments

- boundaries:

  Integer vector: initial boundary positions.

- gene_ranges:

  Data frame with merged gene ranges (columns 'start', 'end').

- chr_end:

  Integer: chromosome end coordinate.

- limit:

  Integer: maximum allowed interval width.

- shift_by:

  Integer: step size (bp) for rightward boundary updates.

- clearance:

  Integer: minimum required clearance (bp) from nearby genes on both
  sides.

- max_iter:

  Integer: maximum number of iterations. Defaults to 100L.

## Value

Integer vector: refined boundary positions after convergence or max
iterations.

## Details

Each iteration:

1.  Shifts all boundaries past overlapping genes

2.  Collects extra boundaries from oversized intervals

3.  Deduplicates and filters boundaries

4.  Shifts again

5.  Checks for convergence (identical to previous iteration)

## Examples

``` r
if (FALSE) { # \dontrun{
  refined <- .refine_boundaries(
    c(100L, 200L), gene_ranges, 300L, 50L, 1L, 1L
  )
} # }
```
