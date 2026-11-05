# Writes a text files with pymol scripts to list significant residues.

Function write a script that can be used in pymol to color structure.
Number of colors and corresponding to them ranges can be defined by
user. Residues are being colored by maximum uptake from significant
peptides per residues.

## Usage

``` r
pymol_script_significant_residue(
  df,
  path = "",
  ranges = c(-Inf, seq(-30, 30, by = 10), Inf),
  alpha = 0.01,
  pvalue = 0.01,
  replicates = 3,
  p_adjust = "none",
  resno_offset = 0,
  pv_cutoff
)
```

## Arguments

- df:

  average data frame. Generated using ave_timepoint() function.

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

pymol script with colors assigned per residues by maximum uptake per
residue

## Examples

``` r
# \donttest{
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- output_tp(file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written
pymol_script_significant_residue(df=a, path=tempdir(), replicates=3, pvalue=0.01,
ranges=c(-Inf,-40, -30,-20,-10, 0,10, 20,30,40, Inf) )

#> NULL
pymol_script_significant_residue(df=a, path=tempdir())

#> NULL
# }
```
