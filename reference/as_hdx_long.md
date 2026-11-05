# Build HDXBoxeR input from a long table

For data from any software: give one row per measurement (peptide,
charge, state, time, replicate) and name the columns that hold each
quantity. The result can be passed to any function that takes
`filepath`.

## Usage

``` r
as_hdx_long(
  data,
  state,
  time,
  replicate,
  start,
  end,
  sequence,
  uptake,
  charge = NULL,
  percent = NULL,
  rt = NULL,
  centroid = NULL,
  time_unit = c("s", "min", "h"),
  undeuterated = 0,
  fd = NULL
)
```

## Arguments

- data:

  data frame with one row per measurement.

- state, time, replicate, start, end, sequence, uptake:

  names of the columns holding the protein state, labeling time,
  replicate (experiment/run name), peptide start and end residue,
  sequence and deuterium uptake in Da.

- charge:

  optional column with the charge state (default: every peptide charge
  1).

- percent:

  optional column with percent deuteration; by default
  `100 * uptake / (length - 1 - prolines after the first residue)`, the
  theoretical maximum HDXExaminer and DynamX use.

- rt:

  optional column with retention time (default 0).

- centroid:

  optional column with the measured centroid m/z (needed only for
  [`backHX_calculations()`](https://mkajano.github.io/HDXBoxeR/reference/backHX_calculations.md)).

- time_unit:

  unit of a numeric `time` column: `"s"` (default), `"min"` or `"h"`.
  Character times are used as they are.

- undeuterated:

  time value of undeuterated controls (default `0`); these rows become
  `"0s"`.

- fd:

  optional time value marking fully deuterated controls (e.g. `"FD"` or
  `Inf`); these rows become `"FD"`.

## Value

data frame of class `hdx_data`.

## See also

[`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)

## Examples

``` r
long <- data.frame(state = rep(c("apo", "holo"), each = 6),
                   t_min = rep(rep(c(1, 10), each = 3), 2),
                   run = rep(paste0("r", 1:3), 4),
                   from = 5, to = 14, seq = "YIAKQRPISF",
                   D = c(1.1, 1.2, 1.15, 2.4, 2.5, 2.45, 0.9, 1.0, 0.95, 2.0, 2.1, 2.05))
d <- as_hdx_long(long, state = "state", time = "t_min", replicate = "run",
                 start = "from", end = "to", sequence = "seq", uptake = "D",
                 time_unit = "min")
output_tp(d)
#> Protein.States used: apo holo
#> Deut.times used: 60.00s 600.00s
#> Number of replicates used: apo=3 holo=3
#> no csv file written
#>   Deut.Time Start End   Sequence Search.RT Charge apo_X..Deut_1 apo_X..Deut_2
#> 1    60.00s     5  14 YIAKQRPISF         0      1           1.1           1.2
#> 2   600.00s     5  14 YIAKQRPISF         0      1           2.4           2.5
#>   apo_X..Deut_3 holo_X..Deut_1 holo_X..Deut_2 holo_X..Deut_3
#> 1          1.15            0.9            1.0           0.95
#> 2          2.45            2.0            2.1           2.05
```
