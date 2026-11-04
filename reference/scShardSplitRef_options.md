# Get or Set scShardSplitRef Options

A convenience function to get or set options used by scShardSplitRef.

## Usage

``` r
scShardSplitRef_options(...)
```

## Arguments

- ...:

  Named options to set, or no arguments to retrieve current values.

## Value

A list of current option values.

## Examples

``` r
# See currently set options for scShardSplitRef
scShardSplitRef_options()
#> $scShardSplitRef.verbosity
#> [1] "verbose"
#> 

# Set verbosity to "quiet" to suppress info/success messages
scShardSplitRef_options(scShardSplitRef.verbosity = "quiet")
scShardSplitRef_options()
#> $scShardSplitRef.verbosity
#> [1] "quiet"
#> 
```
