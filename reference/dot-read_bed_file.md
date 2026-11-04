# Read BED file with validation

Reads a BED-like file and ensures it has the required minimum number of
columns. The first three columns are named 'chr', 'start', 'end'.

## Usage

``` r
.read_bed_file(path, min_cols = 3L)
```

## Arguments

- path:

  Character: file path to BED file.

- min_cols:

  Integer: minimum required number of columns. Defaults to 3L.

## Value

Data frame with columns 'chr', 'start', 'end' plus any additional
columns from the input file.

## Details

Raises an error if the file does not exist or has fewer than `min_cols`
columns. Assumes TAB-delimited format with no header.

## Examples

``` r
if (FALSE) { # \dontrun{
  bed_df <- .read_bed_file("regions.bed", min_cols = 3L)
} # }
```
