# Shift boundaries past genes

Adjusts a vector of boundaries that fall too close to gene intervals by
shifting them rightward. Each boundary is independently shifted until it
satisfies the requested clearance around genes. Results are deduplicated
and filtered.

## Usage

``` r
.shift_boundaries_past_genes(
  boundaries,
  gene_ranges,
  chr_end,
  shift_by = 1L,
  clearance = 1L
)
```

## Arguments

- boundaries:

  Integer vector: boundary positions to adjust.

- gene_ranges:

  Data frame with gene ranges (columns 'start', 'end').

- chr_end:

  Integer: chromosome end coordinate (upper bound).

- shift_by:

  Integer: step size (bp) for rightward boundary updates.

- clearance:

  Integer: minimum required clearance (bp) from nearby genes on both
  sides. Defaults to 1L.

## Value

Integer vector: adjusted, deduplicated, and filtered boundaries. Values
outside (0, chr_end) are removed.

## Examples

``` r
if (FALSE) { # \dontrun{
  shifted <- .shift_boundaries_past_genes(
    c(60L, 150L), gene_ranges, 300L, shift_by = 1L, clearance = 1L
  )
} # }
```
