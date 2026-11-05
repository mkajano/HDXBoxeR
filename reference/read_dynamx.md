# Read a DynamX cluster or state data export

Converts a Waters DynamX export to the HDXBoxeR input table.

## Usage

``` r
read_dynamx(
  path,
  type = c("auto", "cluster", "state"),
  fd = NULL,
  protein = NULL,
  combine_charges = FALSE
)
```

## Arguments

- path:

  path to the csv file.

- type:

  `"auto"` (default), `"cluster"` or `"state"`.

- fd:

  how the fully deuterated control is identified: `NULL` (default,
  none), the name of a DynamX State holding the control (its
  measurements become the `"FD"` control of every other state), or
  `"last"` to use each state's longest exposure as its control (that
  exposure is then not kept as a time point).

- protein:

  which DynamX `Protein` to keep when the file holds more than one.
  Default: all (an error if there are several).

- combine_charges:

  `FALSE` (default) keeps every charge state as its own peptide, as for
  HDXExaminer input. `TRUE` combines charge states of a peptide into one
  intensity-weighted mass per raw file (as in the DynamX state export).

## Value

data frame of class `hdx_data`.

## Details

**Cluster data** (one row per peptide, state, exposure, raw file and
charge; DynamX 2.x or 3.x): each raw `File` is a replicate. `Center` is
the centroid m/z, so the neutral mass is `z * (Center - proton mass)`.
Uptake of a replicate is its mass minus the mean mass of the
undeuterated (`Exposure` 0) measurements of the same peptide, charge and
state. Percent deuteration is `100 * uptake / MaxUptake`.

**State data** (one averaged row per peptide, state and exposure):
`Uptake` is used directly as a single measurement per time point. It
contains no replicates, so functions that need replicate SDs or t-tests
cannot use it; a warning says so.

Exposure is in minutes in DynamX and is converted to seconds (time
labels such as `"600.00s"`; times of 1 s or more are rounded to whole
seconds, so 0.167 min becomes `"10.00s"`). Very short "in-exchange"
controls (e.g. Exposure 0.001 min) are kept as a time point; leave them
out with the `times` argument of
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).
Note that DynamX's `MaxUptake` counts exchangeable amides as length -
1 - prolines, while HDXBoxeR's
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
uses length - 2 - prolines.

## See also

[`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
