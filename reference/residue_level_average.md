# Per-residue uptake as a length-weighted mean across covering peptides

A fast, simple companion to
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md):
for each residue, at each Protein.State/replicate/timepoint, computes a
weighted mean of percent deuteration across every peptide that covers
that residue, with each peptide weighted by
`1 / nb_exch_deut(peptide)` - i.e. shorter, more localized peptides
contribute more than long ones. No kinetic model and no smoothing
between residues are used.

## Usage

``` r
residue_level_average(df, protein_length)
```

## Arguments

- df:

  data.frame, output of
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  called with `percent = TRUE`. Percent deuteration is required here
  (not uptake in Daltons): it is already normalized per peptide length,
  so a weighted average across peptides of different lengths is
  meaningful. Raw uptake is an extensive quantity (scales with peptide
  length) and averaging it directly across peptides of different lengths
  would not be - use
  [`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)
  for uptake-based, length-aware reconstruction.

- protein_length:

  integer, total number of residues in the protein.

## Value

data.frame with one row per residue and `Deut.Time`, and one column per
Protein.State/replicate, giving the weighted mean percent deuteration
across all peptides covering that residue. `NA` where no peptide covers
a residue.

## Details

Unlike
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md),
this does not separate the contribution of overlapping peptides, so a
residue near the edge of a long peptide still carries some signal from
its neighbors. Use it as a quick baseline to compare with
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md).

## See also

[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md),
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md),
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)

## Examples

``` r
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a <- output_tp(file_nm, percent = TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
res <- residue_level_average(a, protein_length = 225)
```
