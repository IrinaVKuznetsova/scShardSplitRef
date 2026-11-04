# Read GTF file with validation

Reads a GTF file and ensures it has the required 9 columns. Comment
lines starting with '#' are automatically skipped.

## Usage

``` r
.read_gtf_file(path)
```

## Arguments

- path:

  Character: file path to GTF file.

## Value

Data frame with standard 9 GTF columns: chr, source, feature, start,
end, score, strand, frame, attribute.

## Details

Raises an error if the file does not exist or has fewer than 9 columns.
Assumes TAB-delimited format with no header.

## Examples

``` r
if (FALSE) { # \dontrun{
  gtf_df <- .read_gtf_file("annotations.gtf")
} # }
```
