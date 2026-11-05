# Provides summary table for all data.sets.

Returns a data frame summarizing each protein state: time points, number
of replicates, number of peptides, sequence coverage, average peptide
length and redundancy.

## Usage

``` r
general_info(filepath, protein_range = NULL)
```

## Arguments

- filepath:

  filepath to the input file. Input file is All_results table from
  HDX_Examiner, where all the fields are marked for export. Can also be
  a DynamX cluster/state export, or a table from
  [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
  /
  [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md).

- protein_range:

  residue range for the coverage, `c(first, last)`. Default: from the
  first to the last residue covered by any peptide in the data (the same
  range for every state).

## Value

Returns summary table.

## Details

Peptides are those found at every time point of the state with at least
two replicates (as in
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md));
each peptide sequence counts once for coverage, length and redundancy
(charge states are not counted twice).

- coverage: percent of residues in `protein_range` covered by at least
  one peptide;

- redundancy: average number of peptides covering a covered residue;

- `# replicates`: replicates per time point (a range when they differ).

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- general_info(file_nm)
general_info(file_nm, protein_range = c(1, 225))
#>         Protein.State                      Timepoints # replicates # peptides
#> Unbound       Unbound 3.00s 60.00s 1800.00s 72000.00s            3         37
#> bound           bound 3.00s 60.00s 1800.00s 72000.00s            3         34
#>         peptide coverage % <Peptide Length> <Redundancy>
#> Unbound              77.78            13.86         2.93
#> bound                77.78            13.29         2.58
```
