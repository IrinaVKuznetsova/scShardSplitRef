# Filter GTF attributes

Extract and keep only specified attributes from GTF attribute column

## Usage

``` r
.filter_attributes(attributes, keep_attributes)
```

## Arguments

- attributes:

  Character vector of GTF attribute strings

- keep_attributes:

  Character vector of attribute names to keep

## Value

Character vector of filtered attribute strings

## Details

Extracts key-value pairs from GTF attribute column and keeps only
specified attributes in the original order.
