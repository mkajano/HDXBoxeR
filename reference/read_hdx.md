# Read HDX-MS data from HDXExaminer or DynamX

Reads an HDX-MS export into the table every HDXBoxeR function works on.
The result can be passed wherever a function asks for `filepath`, e.g.
`output_tp(read_hdx("cluster.csv", fd = "Full Deuteration control"))`.
Passing a file path directly also works; reading once is faster when
several functions use the same data, and is needed for DynamX options
such as `fd`.

## Usage

``` r
read_hdx(
  path,
  format = c("auto", "hdxexaminer", "dynamx_cluster", "dynamx_state"),
  ...
)
```

## Arguments

- path:

  path to the csv file.

- format:

  `"auto"` (default) detects the format from the header; or one of
  `"hdxexaminer"`, `"dynamx_cluster"`, `"dynamx_state"`.

- ...:

  passed to
  [`read_dynamx()`](https://mkajano.github.io/HDXBoxeR/reference/read_dynamx.md)
  (e.g. `fd`, `protein`, `combine_charges`).

## Value

data frame of class `hdx_data` in HDXExaminer's column layout.

## Details

Supported formats (detected from the header line):

- HDXExaminer:

  "All results" table export
  ([`read_hdxexaminer()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdxexaminer.md)).

- DynamX cluster data:

  one row per peptide, state, exposure, raw file and charge
  ([`read_dynamx()`](https://mkajano.github.io/HDXBoxeR/reference/read_dynamx.md));
  DynamX 2.x and 3.x. Raw files are the replicates.

- DynamX state data:

  one averaged row per peptide, state and exposure. It holds no
  replicates, so it can be plotted but the replicate-based significance
  tests need cluster data.

Any other software: build the table from a long data frame with
[`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md).

## See also

[`read_dynamx()`](https://mkajano.github.io/HDXBoxeR/reference/read_dynamx.md),
[`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
d <- read_hdx(file_nm)
a <- output_tp(d)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
```
