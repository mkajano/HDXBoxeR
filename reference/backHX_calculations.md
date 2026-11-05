# Summary of backexchange summary

Returns average and ranges of backexchange. Function calculates as: 1-
(m100%-m0%)/N/Dfact. m0% is the non-deuterated peptide centroid mass,
m100% is the maximally labeled peptide centroid mass, N is the
theoretical number of backbone amides in the peptide and Dfrac is the
fraction of D/H in the labeling buffer used. Function requires
undeuterated and Fully deuterated sets marked in Deut.time as 0s and FD
respectively.

## Usage

``` r
backHX_calculations(filepath, Dfact = 0.85)
```

## Arguments

- filepath:

  filepath to the input file. Input file is All_results table from
  HDX_Examiner, where all the fields are marked for export. Can also be
  a DynamX cluster/state export, or a table from
  [`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
  /
  [`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md).

- Dfact:

  is the fraction of D/H in the labeling buffer used. Default set up to
  0.85

## Value

Returns summary table for backexchange.

## Details

Masses are taken from the measured centroids: HDX Examiner's `Exp Cent`
column is an m/z value, so each centroid is multiplied by the peptide's
`Charge` to convert it to a mass, i.e.
`m100% - m0% = (Exp Cent[FD] - Exp Cent[0s]) * Charge` (in Da). Earlier
versions used the theoretical centroid (`Theor Cent`) without the charge
factor, which overstated back-exchange roughly in proportion to the
charge state. The returned table carries this as its `"method"`
attribute.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a<- backHX_calculations(filepath=file_nm, Dfact=0.85)
```
