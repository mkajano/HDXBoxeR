# Writes a text files with pymol scripts to list significant peptides

Function write a script that can be used in pymol to color structure.
Number of colors and corresponding to them ranges can be defined by
user.

## Usage

``` r
pymol_script_significant_peptide_proc(
  input_proc,
  input_up,
  path = "",
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  alpha = 0.01,
  pvalue = 0.01,
  replicates = 3,
  order.pep = TRUE,
  p_adjust = "none",
  resno_offset = 0,
  pv_cutoff
)
```

## Arguments

- input_proc:

  Dataframe with organized procent deuteration data. Input generated
  using output_tp(, percent=T) function.

- input_up:

  Dataframe with organized deuteration uptake. Input generated using
  output_tp() function.

- path:

  location where the Pymol scripts will be saved. Default `""`: the
  working directory.

- ranges:

  ranges for coloring scheme. Default set to c(-Inf, seq(-30, 30,
  by=10), Inf)

- alpha:

  cutoff for critical interval

- pvalue:

  p-value cutoff here set up to 0.01

- replicates:

  number of replicates in sample. Default set to 3.

- order.pep:

  flag allowing to either order peptide acccording to the peptide length
  (default), or to position in the protein sequence.

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

- resno_offset:

  integer added to every residue number written to the PyMOL script, so
  the data's numbering matches the structure's (e.g. `-20` when the
  structure numbering starts 20 residues later, as with an N-terminal
  tag). Default `0`.

- pv_cutoff:

  deprecated, use `pvalue`.

## Value

pymol script with colors assigned per peptide

## Details

Percent-deuteration counterpart of
[`pymol_script_significant_peptide()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_peptide.md):
significance is tested on the uptake table `input_up` (Da, the
package-wide rule), while peptides are colored by their difference in
percent deuteration, first state minus each state (from `input_proc`);
[`pymol_script_significant_peptide()`](https://mkajano.github.io/HDXBoxeR/reference/pymol_script_significant_peptide.md)
uses the relative change in uptake.

## Examples

``` r
# \donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a_up<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
a_proc<- output_tp(file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
pymol_script_significant_peptide_proc(input_proc=a_proc,
input_up=a_up,  path=tempdir(),replicates=3, alpha=0.01, pvalue=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf), order.pep=TRUE)
#> peptides ordered according to peptide length
#> peptides ordered according to peptide length
#> peptides ordered according to peptide length
#> peptides ordered according to peptide length

#> NULL
# }
```
