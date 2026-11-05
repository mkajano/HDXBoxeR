# Makes input for Extreme for bimodal analysis.

Makes input for Extreme for bimodal analysis.

## Usage

``` r
extreme_input_gap(hm_dir, replicates, timepoints, output_path = "NA")
```

## Arguments

- hm_dir:

  directory in which all the folders which needs to be processed are

- replicates:

  number of replicates in sample

- timepoints:

  lists timepoints used in experiments.

- output_path:

  directory where the output files will be saved, hm_dir default

## Value

Inputs for extreme for all data prepared.

## Examples

``` r
# \donttest{
path_to_folders<-system.file("extdata",  package = "HDXBoxeR")

extreme_input_gap(hm_dir =path_to_folders, replicates = 3,
timepoints =c(3, 60, 1800, 72000), output_path=tempdir())
#> processing:.
#> processing:./0001-0014-MQIFVKTLTGKTIT
# }
```
