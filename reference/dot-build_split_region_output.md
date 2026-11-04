# Build final split-region output data frame

Build final split-region output data frame

## Usage

``` r
.build_split_region_output(
  regions,
  merged_gene_ranges,
  limit,
  shift_by,
  clearance
)
```

## Arguments

- regions:

  BED-like regions data frame.

- merged_gene_ranges:

  Named list of merged gene ranges by chromosome.

- limit:

  Maximum interval width.

- shift_by:

  Boundary shift step size.

- clearance:

  Boundary clearance distance.

## Value

Data frame of final split intervals.
