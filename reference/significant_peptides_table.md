# Interactive table of peptides, linked to the volcano/fissure_vent/Woods plots

[`DT::datatable()`](https://rdrr.io/pkg/DT/man/datatable.html) of every
peptide in the comparison (sortable/searchable), with a `significant`
column from the same hybrid criterion as the other interactive plots.
Selecting rows highlights the corresponding peptides on any other plot
built from the same `shared_data`.

## Usage

``` r
significant_peptides_table(
  df,
  state1,
  state2,
  times,
  replicates = 3,
  alpha = 0.01,
  pvalue = 0.01,
  shared_data = NULL,
  p_adjust = "none"
)
```

## Arguments

- df:

  output of
  [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md).

- state1:

  First protein state, matched as in
  [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md).
  Defaults to the first state found in `df`.

- state2:

  Second protein state. Defaults to the second state found.

- times:

  a single Deut.Time value to plot.

- replicates:

  number of replicates. Default `3`.

- alpha:

  alpha for the pooled critical interval. Default `0.01`. May be a
  vector (paired element-wise with `pvalue`, recycled to a common
  length) - with more than one combination, the plot gets a dropdown
  (native `plotly` `updatemenus`, no R session needed) that swaps
  between the pre-computed colorings/thresholds for each
  `alpha`/`pvalue` pair. Peptide *positions* never move (delta and
  p-value don't depend on either parameter) - only which points/lines
  count as significant.

- pvalue:

  per-peptide significance cutoff. Default `0.01`. See `alpha` for
  vector/dropdown behavior.

- shared_data:

  optional
  [`crosstalk::SharedData`](https://rdrr.io/pkg/crosstalk/man/SharedData.html)
  object (as built internally by
  [`interactive_fissure_explorer()`](https://mkajano.github.io/HDXBoxeR/reference/interactive_fissure_explorer.md))
  so this plot's selection is linked to other plots sharing the same
  object. If `NULL` (default), a new one is built and this plot works
  standalone. The linked *key* (which peptide is which) always follows
  `shared_data`; only the significance coloring/dropdown below is local
  to this plot.

- p_adjust:

  multiple-testing correction of the per-peptide Welch p-values, passed
  to [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html):
  `"none"` (default) or e.g. `"BH"`, applied across the peptides of the
  comparison.

## Value

a [`DT::datatable`](https://rdrr.io/pkg/DT/man/datatable.html)
htmlwidget.

## Examples

``` r
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)
    && requireNamespace("DT", quietly = TRUE)) {
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a <- output_tp(file_nm)
significant_peptides_table(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
}
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written

{"x":{"crosstalkOptions":{"key":["47_63_VVPGPAEHPLQYNYTFW_z2","89_97_ASVEQFWRF_z2","114_124_HLFKEGIKPMW_z2","136_149_IIRLRKGLASRCWE_z3","116_135_FKEGIKPMWEDDANKNGGKW_z3","62_77_FWYSRRTPGRPTSSQS_z3","194_200_ARIRDTL_z2","82_88_IKQIGTF_z2","116_124_FKEGIKPMW_z2","169_177_AVVSVRFQE_z2","95_115_WRFYSHMVRPGDLTGHSDFHL_z3","63_77_WYSRRTPGRPTSSQS_z3","98_115_YSHMVRPGDLTGHSDFHL_z2","154_161_AMLGEQFM_z2","213_225_EYKTHTDSIKMPG_z2","78_88_YEQNIKQIGTF_z2","213_220_EYKTHTDS_z2","114_135_HLFKEGIKPMWEDDANKNGGKW_z3","136_147_IIRLRKGLASRC_z3","40_58_GPLGSKAVVPGPAEHPLQY_z2","95_107_WRFYSHMVRPGDL_z3","40_59_GPLGSKAVVPGPAEHPLQYN_z2","183_193_WNKTASDQATT_z2","214_220_YKTHTDS_z2","64_77_YSRRTPGRPTSSQS_z2","201_213_RRVLNLPPNTIME_z2","183_200_WNKTASDQATTARIRDTL_z2","176_183_QEDIISIW_z2","178_193_DIISIWNKTASDQATT_z2","140_147_RKGLASRC_z2","176_193_QEDIISIWNKTASDQATT_z2","182_193_IWNKTASDQATT_z2"],"group":"hdx_table_nbesrxjl"},"filter":"top","vertical":false,"filterHTML":"<tr>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"character\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"40\" data-max=\"214\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"integer\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"58\" data-max=\"225\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"-1.328\" data-max=\"0.254\" data-scale=\"3\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"number\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"display: none;position: absolute;width: 200px;opacity: 1\">\n      <div data-min=\"3.57e-07\" data-max=\"0.464\" data-scale=\"9\"><\/div>\n      <span style=\"float: left;\"><\/span>\n      <span style=\"float: right;\"><\/span>\n    <\/div>\n  <\/td>\n  <td data-type=\"logical\" style=\"vertical-align: top;\">\n    <div class=\"form-group has-feedback\" style=\"margin-bottom: auto;\">\n      <input type=\"search\" placeholder=\"All\" class=\"form-control\" style=\"width: 100%;\"/>\n      <span class=\"glyphicon glyphicon-remove-circle form-control-feedback\"><\/span>\n    <\/div>\n    <div style=\"width: 100%; display: none;\">\n      <select multiple=\"multiple\" style=\"width: 100%;\" data-options=\"[&quot;true&quot;,&quot;false&quot;]\"><\/select>\n    <\/div>\n  <\/td>\n<\/tr>","data":[["47_63_VVPGPAEHPLQYNYTFW_z2","89_97_ASVEQFWRF_z2","114_124_HLFKEGIKPMW_z2","136_149_IIRLRKGLASRCWE_z3","116_135_FKEGIKPMWEDDANKNGGKW_z3","62_77_FWYSRRTPGRPTSSQS_z3","194_200_ARIRDTL_z2","82_88_IKQIGTF_z2","116_124_FKEGIKPMW_z2","169_177_AVVSVRFQE_z2","95_115_WRFYSHMVRPGDLTGHSDFHL_z3","63_77_WYSRRTPGRPTSSQS_z3","98_115_YSHMVRPGDLTGHSDFHL_z2","154_161_AMLGEQFM_z2","213_225_EYKTHTDSIKMPG_z2","78_88_YEQNIKQIGTF_z2","213_220_EYKTHTDS_z2","114_135_HLFKEGIKPMWEDDANKNGGKW_z3","136_147_IIRLRKGLASRC_z3","40_58_GPLGSKAVVPGPAEHPLQY_z2","95_107_WRFYSHMVRPGDL_z3","40_59_GPLGSKAVVPGPAEHPLQYN_z2","183_193_WNKTASDQATT_z2","214_220_YKTHTDS_z2","64_77_YSRRTPGRPTSSQS_z2","201_213_RRVLNLPPNTIME_z2","183_200_WNKTASDQATTARIRDTL_z2","176_183_QEDIISIW_z2","178_193_DIISIWNKTASDQATT_z2","140_147_RKGLASRC_z2","176_193_QEDIISIWNKTASDQATT_z2","182_193_IWNKTASDQATT_z2"],["VVPGPAEHPLQYNYTFW","ASVEQFWRF","HLFKEGIKPMW","IIRLRKGLASRCWE","FKEGIKPMWEDDANKNGGKW","FWYSRRTPGRPTSSQS","ARIRDTL","IKQIGTF","FKEGIKPMW","AVVSVRFQE","WRFYSHMVRPGDLTGHSDFHL","WYSRRTPGRPTSSQS","YSHMVRPGDLTGHSDFHL","AMLGEQFM","EYKTHTDSIKMPG","YEQNIKQIGTF","EYKTHTDS","HLFKEGIKPMWEDDANKNGGKW","IIRLRKGLASRC","GPLGSKAVVPGPAEHPLQY","WRFYSHMVRPGDL","GPLGSKAVVPGPAEHPLQYN","WNKTASDQATT","YKTHTDS","YSRRTPGRPTSSQS","RRVLNLPPNTIME","WNKTASDQATTARIRDTL","QEDIISIW","DIISIWNKTASDQATT","RKGLASRC","QEDIISIWNKTASDQATT","IWNKTASDQATT"],[47,89,114,136,116,62,194,82,116,169,95,63,98,154,213,78,213,114,136,40,95,40,183,214,64,201,183,176,178,140,176,182],[63,97,124,149,135,77,200,88,124,177,115,77,115,161,225,88,220,135,147,58,107,59,193,220,77,213,200,183,193,147,193,193],[-1.328,-1.129,-0.364,-1.007,-0.731,-0.483,-0.394,-0.275,-0.38,-0.357,-0.711,0.254,-0.495,-0.148,-0.239,-0.114,-0.308,-0.492,-0.212,-0.141,-0.157,0.129,0.124,0.08400000000000001,-0.15,-0.066,-0.346,-0.07099999999999999,-0.057,-0.031,-0.091,0],[3.57e-07,6.600000000000001e-05,8.610000000000001e-05,0.000102,0.000313,0.000409,0.000556,0.000753,0.00131,0.00143,0.0021,0.00273,0.00309,0.00315,0.00349,0.008330000000000001,0.0108,0.0152,0.0166,0.037,0.0435,0.0586,0.0761,0.09080000000000001,0.102,0.117,0.182,0.222,0.339,0.351,0.464,null],[true,true,true,true,true,true,true,true,true,true,true,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th>peptide_id<\/th>\n      <th>Sequence<\/th>\n      <th>Start<\/th>\n      <th>End<\/th>\n      <th>Delta<\/th>\n      <th>p_value<\/th>\n      <th>Significant<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"pageLength":10,"columnDefs":[{"visible":false,"targets":0},{"className":"dt-right","targets":[2,3,4,5]},{"name":"peptide_id","targets":0},{"name":"Sequence","targets":1},{"name":"Start","targets":2},{"name":"End","targets":3},{"name":"Delta","targets":4},{"name":"p_value","targets":5},{"name":"Significant","targets":6}],"order":[],"autoWidth":false,"orderClasses":false,"orderCellsTop":true},"selection":{"mode":"multiple","selected":null,"target":"row","selectable":null}},"evals":[],"jsHooks":[]}
```
