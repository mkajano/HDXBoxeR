# Per-residue deconvolution of HDX-MS uptake from overlapping peptides

Estimates a single-residue-resolution deuterium uptake profile from
overlapping peptide-level data, by solving a constrained least-squares
problem rather than simply averaging peptides that cover a residue. Each
peptide's observed uptake (in Daltons) is modelled as the sum of the
per-residue uptake of every exchangeable residue it covers (the first
two residues of each peptide and all non-first Prolines are excluded,
matching the convention used in
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)).
The resulting system is solved per Protein State, replicate, and
timepoint with residue-level uptake constrained to `\[0, 1\]` (a residue
cannot be less than fully protected or more than fully deuterated). With
typical overlapping-peptide coverage (redundancy \> ~5x), this recovers
sharper, more localized residue-level resolution than a plain
per-residue mean across covering peptides, particularly at region
boundaries where long peptides would otherwise blend signal from
neighboring residues into a region's apparent value.

## Usage

``` r
residue_deconvolution(df, protein_length, lambda_smooth = 0.1)
```

## Arguments

- df:

  data.frame, output of
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  called with `percent = FALSE`. Uptake in Daltons is required because
  it is additive across residues; percent deuteration is not (it depends
  on each peptide's own number of exchangeable residues) and cannot be
  used directly for deconvolution.

- protein_length:

  integer, total number of residues in the protein. Used to set the
  residue index range for the deconvolution.

- lambda_smooth:

  non-negative numeric, default `0.1`. Strength of an L2 penalty on the
  difference between adjacent covered residues (only applied between
  residues that are actually consecutive in sequence - a gap in coverage
  breaks the chain). `0` reproduces the original, unregularized
  per-timepoint fit. Larger values pull neighboring residues toward the
  same value unless the peptide data specifically distinguish them. The
  default of `0.1` is a modest starting point, not a universally correct
  value - the right strength depends on peptide redundancy and protein
  size. On the example data, values between 0 and 1 predicted left-out
  peptides equally well.

## Value

data.frame with one row per residue and `Deut.Time`, and one column per
Protein.State/replicate (same naming as the `X..Deut` columns in `df`),
giving the deconvolved per-residue fractional uptake (0-1). Residues not
covered by any peptide at a given timepoint are `NA`.

## What it does not do

Each time point is solved on its own. There is no kinetic model and no
link between time points. `lambda_smooth > 0` adds a penalty on the
difference between neighboring residues within one time point.

## See also

[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md),
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md),
[`coverage_residue()`](https://mkajano.github.io/HDXBoxeR/reference/coverage_residue.md),
[`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a <- output_tp(file_nm, percent = FALSE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
res <- residue_deconvolution(a, protein_length = 225)
res_unsmoothed <- residue_deconvolution(a, protein_length = 225, lambda_smooth = 0)
```
