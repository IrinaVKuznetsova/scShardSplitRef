# Write BED file

Writes a data frame to a TAB-delimited BED format file without quotes or
row names.

## Usage

``` r
.write_bed_file(df, path)
```

## Arguments

- df:

  Data frame with 'chr', 'start', 'end' columns.

- path:

  Character: output file path.

## Value

Called for its side effects, writes a file to local disk.

## Examples

``` r
if (FALSE) { # \dontrun{
  .write_bed_file(split_regions_df, "output.bed")
} # }
```
