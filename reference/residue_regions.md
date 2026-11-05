# Cluster significant peptides into contiguous destabilized/protected regions

Takes a peptide-level comparison between two states at one timepoint and
derives contiguous residue regions that are consistently destabilized
(more exposed) or protected (less exposed), rather than reporting a flat
list of significant peptides. A peptide counts as significant by the
same rule as everywhere else in the package: Welch p-value below
`pvalue` (raw, or adjusted with `p_adjust`) AND absolute difference
above the pooled Hageman & Weis critical interval. A residue enters a
region only once at least `min_peptides` distinct significant peptides
of the *same direction* cover it, so a single peptide alone cannot
create a region. Charge states of the same peptide count once.

## Usage

``` r
residue_regions(
  df,
  state1,
  state2,
  times,
  replicates = 3,
  pvalue = 0.01,
  p_adjust = "none",
  min_peptides = 2,
  gap_tolerance = 2,
  alpha = 0.01,
  fdr_cutoff
)
```

## Arguments

- df:

  output of
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
  called with `percent = FALSE`, i.e. deuteration uptake in Daltons -
  this function works in Da so that `delta` is additive across residues,
  same requirement as
  [`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md).

- state1:

  First protein state, matched the same way as in
  [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md).

- state2:

  Second protein state, same matching rule.

- times:

  a single Deut.Time value to derive regions for (e.g. `"1800.00s"`).
  Call this function once per timepoint of interest and
  [`rbind()`](https://rdrr.io/r/base/cbind.html) the results if you want
  regions across multiple timepoints.

- replicates:

  number of replicates in sample, only used for tables without named
  replicate columns. Default set to 3.

- pvalue:

  p-value cutoff for a peptide to count as significant. Default `0.01`.

- p_adjust:

  multiple-testing correction of the per-peptide Welch p-values, passed
  to [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html):
  `"none"` (default) or e.g. `"BH"`, applied across the peptides of the
  comparison.

- min_peptides:

  minimum number of distinct, same-direction, significant peptides that
  must cover a residue before that residue is included in a region.
  Default `2`.

- gap_tolerance:

  maximum number of consecutive non-covered residues allowed inside one
  region before it is split into two. Default `2`.

- alpha:

  alpha of the pooled Hageman & Weis critical interval
  ([`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md)),
  the same threshold
  [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
  draws as its grey band. Default `0.01`.

- fdr_cutoff:

  deprecated. Earlier versions called a peptide significant when its
  BH-adjusted p-value was below `fdr_cutoff`, without the
  critical-interval condition. Supplying it reproduces that rule, with a
  warning.

## Value

data.frame with one row per region: `direction`
("destabilized"/"protected"), `region_start`, `region_end`,
`n_residues`, `n_peptides` (distinct significant, same-direction
peptides overlapping the region), `mean_per_residue_delta` (Da per
exchangeable amide), and `pooled_CI_Da` (the pooled critical interval
for this comparison, repeated on every row for reference).

## Details

The per-residue effect reported for a region is each covering peptide's
difference divided by its number of exchangeable amides
([`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)),
averaged over the peptides covering each residue.

## See also

[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md),
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md),
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
residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#>   direction region_start region_end n_residues n_peptides
#> 1 protected           62         63          2          2
#> 2 protected           95        124         30          6
#>   mean_per_residue_delta pooled_CI_Da
#> 1            -0.07548611    0.2392609
#> 2            -0.04715929    0.2392609
```
