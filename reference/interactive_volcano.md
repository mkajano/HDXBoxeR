# Interactive volcano plot, linked to the fissure_vent/Woods/table views

Plotly version of
[`plots_vol_tp()`](https://mkajano.github.io/HDXBoxeR/reference/plots_vol_tp.md)'s
delta-vs-p-value display: x = state2

- state1 (Da or the units of `df`), y = -log10(p-value), with the same
  hybrid criterion (Hageman & Weis pooled CI + per-peptide Welch
  p-value) used everywhere else in this package. Hover shows the peptide
  sequence, range, both states' values, delta, and p-value.

## Usage

``` r
interactive_volcano(
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

a `plotly` htmlwidget.

## Examples

``` r
if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("crosstalk", quietly = TRUE)) {
file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
a <- output_tp(file_nm)
interactive_volcano(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
interactive_volcano(a, state1 = "Unbound", state2 = "bound", times = "1800.00s",
                     alpha = c(0.01, 0.05), pvalue = c(0.01, 0.05))
}
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 29 peptide row(s) dropped: 29 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: Unbound=3 bound=3
#> no csv file written

{"x":{"visdat":{"1b024b8ef5d5":["function () ","plotlyVisDat"]},"cur_data":"1b024b8ef5d5","attrs":{"1b024b8ef5d5":{"x":{},"y":{},"text":{},"hoverinfo":"text","mode":"markers","marker":{"size":9,"opacity":0.75,"color":["grey60","grey60","firebrick","firebrick","firebrick","grey60","grey60","firebrick","firebrick","grey60","firebrick","firebrick","firebrick","grey60","firebrick","firebrick","grey60","firebrick","grey60","grey60","firebrick","grey60","grey60","grey60","grey60","grey60","firebrick","grey60","grey60","firebrick","grey60"]},"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"xaxis":{"domain":[0,1],"automargin":true,"title":"delta (state2 - state1)"},"yaxis":{"domain":[0,1],"automargin":true,"title":"-log10(p-value)"},"shapes":[{"type":"line","x0":0,"x1":0,"xref":"x","y0":-1,"y1":1,"yref":"paper","line":{"color":"grey80","dash":"dot"}},{"type":"line","x0":-1,"x1":1,"xref":"paper","y0":2,"y1":2,"yref":"y","line":{"color":"grey60","dash":"dash"}}],"updatemenus":[{"type":"dropdown","active":0,"x":1.02,"y":1,"xanchor":"left","pad":{"r":10,"t":10},"buttons":[{"method":"update","label":"alpha=0.01, p=0.01","args":[{"marker.color":[["grey60","grey60","firebrick","firebrick","firebrick","grey60","grey60","firebrick","firebrick","grey60","firebrick","firebrick","firebrick","grey60","firebrick","firebrick","grey60","firebrick","grey60","grey60","firebrick","grey60","grey60","grey60","grey60","grey60","firebrick","grey60","grey60","firebrick","grey60"]]},{"shapes":[{"type":"line","x0":0,"x1":0,"xref":"x","y0":-1,"y1":1,"yref":"paper","line":{"color":"grey80","dash":"dot"}},{"type":"line","x0":-1,"x1":1,"xref":"paper","y0":2,"y1":2,"yref":"y","line":{"color":"grey60","dash":"dash"}}]}]},{"method":"update","label":"alpha=0.05, p=0.05","args":[{"marker.color":[["grey60","grey60","firebrick","firebrick","firebrick","grey60","grey60","firebrick","firebrick","firebrick","firebrick","firebrick","firebrick","firebrick","firebrick","firebrick","firebrick","firebrick","grey60","firebrick","firebrick","grey60","grey60","grey60","grey60","grey60","firebrick","grey60","firebrick","firebrick","grey60"]]},{"shapes":[{"type":"line","x0":0,"x1":0,"xref":"x","y0":-1,"y1":1,"yref":"paper","line":{"color":"grey80","dash":"dot"}},{"type":"line","x0":-1,"x1":1,"xref":"paper","y0":1.3010299956639813,"y1":1.3010299956639813,"yref":"y","line":{"color":"grey60","dash":"dash"}}]}]}]}],"dragmode":"zoom","hovermode":"closest","showlegend":false},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"highlight":{"on":"plotly_hover","off":"plotly_doubleclick","persistent":false,"dynamic":false,"color":null,"selectize":false,"defaultValues":null,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0,"ctGroups":["hdx_volcano_vjubqtaf"]},"data":[{"x":[-0.14100000000000001,0.1293333333333333,-1.3283333333333331,-0.4833333333333325,0.25399999999999956,-0.14966666666666661,-0.11366666666666614,-0.27499999999999947,-1.1289999999999996,-0.15733333333333377,-0.71066666666666656,-0.49500000000000011,-0.36366666666666614,-0.49233333333333462,-0.3796666666666666,-0.73066666666666791,-0.21200000000000063,-1.0069999999999997,-0.030666666666666842,-0.14833333333333387,-0.35666666666666647,-0.071333333333333332,-0.091000000000000192,-0.057333333333333236,0.12366666666666593,-0.34600000000000009,-0.39366666666666683,-0.065999999999999837,-0.30833333333333313,-0.23933333333333362,0.084333333333333371],"y":[1.4312134261800504,1.2317840640167312,6.4478609387401891,3.3878178673061896,2.5637870202904738,0.99083802688732792,2.0794882176784304,3.1233705728279828,4.1801686199759791,1.3611209788319736,2.6776053047084796,2.5094149008104436,4.0649106122406202,1.8181636705214181,2.8841052615485143,3.5043180761706694,1.7793798516626991,3.99163268840418,0.45510492799829932,2.5014092993404251,2.8432656623411265,0.6542785344535893,0.33393213833625873,0.46992020774452931,1.1188498983340698,0.7392306466609404,3.2547046191304174,0.93336011101893257,1.9669972079922653,2.4568614363398376,1.0416967138879638],"text":["GPLGSKAVVPGPAEHPLQY<br />Start-End: 40-58<br />Unbound: 9.32  bound: 9.18<br />delta: -0.14   p = 0.037","GPLGSKAVVPGPAEHPLQYN<br />Start-End: 40-59<br />Unbound: 9.90  bound: 10.03<br />delta: 0.13   p = 0.0586","VVPGPAEHPLQYNYTFW<br />Start-End: 47-63<br />Unbound: 6.96  bound: 5.63<br />delta: -1.33   p = 3.57e-07  [significant]","FWYSRRTPGRPTSSQS<br />Start-End: 62-77<br />Unbound: 9.30  bound: 8.82<br />delta: -0.48   p = 0.000409  [significant]","WYSRRTPGRPTSSQS<br />Start-End: 63-77<br />Unbound: 8.23  bound: 8.48<br />delta: 0.25   p = 0.00273  [significant]","YSRRTPGRPTSSQS<br />Start-End: 64-77<br />Unbound: 7.78  bound: 7.63<br />delta: -0.15   p = 0.102","YEQNIKQIGTF<br />Start-End: 78-88<br />Unbound: 6.70  bound: 6.59<br />delta: -0.11   p = 0.00833","IKQIGTF<br />Start-End: 82-88<br />Unbound: 4.06  bound: 3.79<br />delta: -0.27   p = 0.000753  [significant]","ASVEQFWRF<br />Start-End: 89-97<br />Unbound: 5.76  bound: 4.63<br />delta: -1.13   p = 6.6e-05  [significant]","WRFYSHMVRPGDL<br />Start-End: 95-107<br />Unbound: 5.95  bound: 5.79<br />delta: -0.16   p = 0.0435","WRFYSHMVRPGDLTGHSDFHL<br />Start-End: 95-115<br />Unbound: 7.93  bound: 7.22<br />delta: -0.71   p = 0.0021  [significant]","YSHMVRPGDLTGHSDFHL<br />Start-End: 98-115<br />Unbound: 5.81  bound: 5.32<br />delta: -0.50   p = 0.00309  [significant]","HLFKEGIKPMW<br />Start-End: 114-124<br />Unbound: 4.31  bound: 3.95<br />delta: -0.36   p = 8.61e-05  [significant]","HLFKEGIKPMWEDDANKNGGKW<br />Start-End: 114-135<br />Unbound: 9.52  bound: 9.03<br />delta: -0.49   p = 0.0152","FKEGIKPMW<br />Start-End: 116-124<br />Unbound: 4.25  bound: 3.87<br />delta: -0.38   p = 0.00131  [significant]","FKEGIKPMWEDDANKNGGKW<br />Start-End: 116-135<br />Unbound: 9.31  bound: 8.58<br />delta: -0.73   p = 0.000313  [significant]","IIRLRKGLASRC<br />Start-End: 136-147<br />Unbound: 7.59  bound: 7.38<br />delta: -0.21   p = 0.0166","IIRLRKGLASRCWE<br />Start-End: 136-149<br />Unbound: 8.76  bound: 7.76<br />delta: -1.01   p = 0.000102  [significant]","RKGLASRC<br />Start-End: 140-147<br />Unbound: 4.83  bound: 4.80<br />delta: -0.03   p = 0.351","AMLGEQFM<br />Start-End: 154-161<br />Unbound: 3.18  bound: 3.03<br />delta: -0.15   p = 0.00315","AVVSVRFQE<br />Start-End: 169-177<br />Unbound: 3.58  bound: 3.23<br />delta: -0.36   p = 0.00143  [significant]","QEDIISIW<br />Start-End: 176-183<br />Unbound: 0.21  bound: 0.14<br />delta: -0.07   p = 0.222","QEDIISIWNKTASDQATT<br />Start-End: 176-193<br />Unbound: 6.92  bound: 6.83<br />delta: -0.09   p = 0.464","DIISIWNKTASDQATT<br />Start-End: 178-193<br />Unbound: 6.78  bound: 6.72<br />delta: -0.06   p = 0.339","WNKTASDQATT<br />Start-End: 183-193<br />Unbound: 6.41  bound: 6.53<br />delta: 0.12   p = 0.0761","WNKTASDQATTARIRDTL<br />Start-End: 183-200<br />Unbound: 10.79  bound: 10.44<br />delta: -0.35   p = 0.182","ARIRDTL<br />Start-End: 194-200<br />Unbound: 3.31  bound: 2.92<br />delta: -0.39   p = 0.000556  [significant]","RRVLNLPPNTIME<br />Start-End: 201-213<br />Unbound: 7.11  bound: 7.05<br />delta: -0.07   p = 0.117","EYKTHTDS<br />Start-End: 213-220<br />Unbound: 2.44  bound: 2.13<br />delta: -0.31   p = 0.0108","EYKTHTDSIKMPG<br />Start-End: 213-225<br />Unbound: 4.83  bound: 4.59<br />delta: -0.24   p = 0.00349  [significant]","YKTHTDS<br />Start-End: 214-220<br />Unbound: 2.04  bound: 2.13<br />delta: 0.08   p = 0.0908"],"hoverinfo":["text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text","text"],"mode":"markers","marker":{"color":["grey60","grey60","firebrick","firebrick","firebrick","grey60","grey60","firebrick","firebrick","grey60","firebrick","firebrick","firebrick","grey60","firebrick","firebrick","grey60","firebrick","grey60","grey60","firebrick","grey60","grey60","grey60","grey60","grey60","firebrick","grey60","grey60","firebrick","grey60"],"size":9,"opacity":0.75,"line":{"color":"rgba(31,119,180,1)"}},"type":"scatter","key":["40_58_GPLGSKAVVPGPAEHPLQY_z2","40_59_GPLGSKAVVPGPAEHPLQYN_z2","47_63_VVPGPAEHPLQYNYTFW_z2","62_77_FWYSRRTPGRPTSSQS_z3","63_77_WYSRRTPGRPTSSQS_z3","64_77_YSRRTPGRPTSSQS_z2","78_88_YEQNIKQIGTF_z2","82_88_IKQIGTF_z2","89_97_ASVEQFWRF_z2","95_107_WRFYSHMVRPGDL_z3","95_115_WRFYSHMVRPGDLTGHSDFHL_z3","98_115_YSHMVRPGDLTGHSDFHL_z2","114_124_HLFKEGIKPMW_z2","114_135_HLFKEGIKPMWEDDANKNGGKW_z3","116_124_FKEGIKPMW_z2","116_135_FKEGIKPMWEDDANKNGGKW_z3","136_147_IIRLRKGLASRC_z3","136_149_IIRLRKGLASRCWE_z3","140_147_RKGLASRC_z2","154_161_AMLGEQFM_z2","169_177_AVVSVRFQE_z2","176_183_QEDIISIW_z2","176_193_QEDIISIWNKTASDQATT_z2","178_193_DIISIWNKTASDQATT_z2","183_193_WNKTASDQATT_z2","183_200_WNKTASDQATTARIRDTL_z2","194_200_ARIRDTL_z2","201_213_RRVLNLPPNTIME_z2","213_220_EYKTHTDS_z2","213_225_EYKTHTDSIKMPG_z2","214_220_YKTHTDS_z2"],"set":"hdx_volcano_vjubqtaf","error_y":{"color":"rgba(31,119,180,1)"},"error_x":{"color":"rgba(31,119,180,1)"},"line":{"color":"rgba(31,119,180,1)"},"xaxis":"x","yaxis":"y","_isNestedKey":false,"frame":null}],"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}
```
