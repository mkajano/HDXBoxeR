# Roll up derived regions against named residue subsets

A companion to
[`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md):
that function always scans every residue covered by any tested peptide -
the unrestricted, data-driven "total" picture, deliberately not
pre-restricted to any subset you might already suspect matters. This
function takes that unrestricted result and a named list of residue
subsets (e.g. structural domains) and reports, per subset, how much of
it the data-driven regions actually cover - a cross-check for a
structural hypothesis against the unconstrained result, not a
replacement for it. Call
[`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md)
once (it does the expensive per-peptide testing and clustering), then
this function as many times as you like against different subset
definitions - it's a cheap rollup over an already-computed result.

## Usage

``` r
residue_regions_subset_summary(region_df, regions)
```

## Arguments

- region_df:

  the data.frame returned by
  [`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md).

- regions:

  named list of `c(start, end)` residue ranges to check against
  `region_df`, e.g. `list(NTR = c(1, 85), ACD = c(86, 165))`.

## Value

data.frame with one row per named subset: `subset`, `n_residues` (subset
size), `n_residues_significant` (residues in the subset covered by any
region in `region_df`), `frac_significant`, `n_regions_overlapping`, and
`directions` (comma-separated set of `direction` values among the
overlapping regions, `""` if none).

## See also

[`residue_regions()`](https://mkajano.github.io/HDXBoxeR/reference/residue_regions.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a <- output_tp(file_nm, percent = FALSE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
rr <- residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
residue_regions_subset_summary(rr, regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#>   subset n_residues n_residues_significant frac_significant
#> 1  nterm        100                      8            0.080
#> 2  cterm        125                     24            0.192
#>   n_regions_overlapping directions
#> 1                     2  protected
#> 2                     1  protected
```
