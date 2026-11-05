# Per-peptide quality checks of the input data

Lists, for every peptide (start, end, sequence, charge) in every protein
state, simple checks of the controls and uptake values. Nothing is
removed from the data; use the table to decide what to exclude.

## Usage

``` r
peptide_qc(filepath)
```

## Arguments

- filepath:

  input file or table, as for
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

## Value

data frame with one row per state and peptide.

## Details

Checks:

- `has_undeut`:

  an undeuterated (`"0s"`) control is present. Without it HDXExaminer
  cannot compute uptake and reports 0.

- `all_zero_uptake`:

  uptake is exactly 0 at every labeling time and in the FD control - in
  practice a missing reference, not a real measurement of no exchange.

- `negative_uptake`:

  at least one replicate has uptake below 0.

- `fd_uptake`:

  mean uptake (Da) of the fully deuterated (`"FD"`) control, `NA` if
  there is none.

- `max_uptake`:

  theoretical maximum: length - 1 - prolines after the first residue (as
  HDXExaminer and DynamX count it).

- `fd_fraction`:

  `fd_uptake / max_uptake`; 1 minus this is the peptide's back-exchange
  (without correction for the D2O fraction).

- `fd_missing`, `fd_nonpositive`, `fd_above_max`:

  no FD value; FD uptake \<= 0; FD uptake above the theoretical maximum.

- `issues`:

  the checks that failed, as text (empty if none).

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
qc <- peptide_qc(file_nm)
qc[qc$issues != "", ]
#>    Protein.State Start End     Sequence Charge has_undeut all_zero_uptake
#> 30       Unbound   182 193 IWNKTASDQATT      2      FALSE            TRUE
#> 64         bound   182 193 IWNKTASDQATT      2      FALSE            TRUE
#>    negative_uptake fd_uptake max_uptake fd_fraction fd_missing fd_nonpositive
#> 30           FALSE         0         11           0      FALSE           TRUE
#> 64           FALSE         0         11           0      FALSE           TRUE
#>    fd_above_max                                                          issues
#> 30        FALSE no undeuterated control; uptake 0 at every time; FD uptake <= 0
#> 64        FALSE no undeuterated control; uptake 0 at every time; FD uptake <= 0
```
