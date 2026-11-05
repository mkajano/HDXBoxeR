# Input Formats, Replicates and Normalization

This vignette shows how to load data, how replicates are handled, and
which value scales are available.

``` r

library(HDXBoxeR)
```

## Reading data

Every function that takes `filepath` also takes a table that is already
loaded. Load the data once with
[`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
and pass it on.

### HDXExaminer

Export from HDXExaminer: Tools \> Export \> All results table (see
[`?All_results_table`](https://mkajano.github.io/HDXBoxeR/reference/All_results_table.md)).

``` r

f3 <- system.file("extdata", "three_states.csv", package = "HDXBoxeR")
d3 <- read_hdx(f3)
table(d3$Protein.State, d3$Deut.Time)
#>              
#>                0s 1800.00s 4.00s 60.00s  FD
#>   D109H_pH7.5 150      200   150    200  85
#>   R120G_ph7.5  99      198   150    150  75
#>   WT_pH7.5    149      200   200    200 132
```

### Waters DynamX

[`read_hdx()`](https://mkajano.github.io/HDXBoxeR/reference/read_hdx.md)
also reads DynamX cluster data and state data (DynamX 2.x and 3.x).

- Cluster data: each raw file is one replicate. Uptake is calculated
  from the centroid mass.
- State data: one average value per peptide, without replicates. It can
  be plotted, but the significance tests need cluster data.
- DynamX exposure times (minutes) become seconds: 0.5 min is `"30.00s"`.

``` r

d <- read_hdx("cluster_data.csv",
              fd = "Full Deuteration control",  # state with the FD control, or "last"
              combine_charges = FALSE)          # TRUE merges charge states
```

### Other software

[`as_hdx_long()`](https://mkajano.github.io/HDXBoxeR/reference/as_hdx_long.md)
takes a table with one row per measurement. Name the columns that hold
each value:

``` r

long <- data.frame(state = rep(c("apo", "holo"), each = 6),
                   t_min = rep(rep(c(1, 10), each = 3), 2),
                   run = rep(paste0("r", 1:3), 4),
                   from = 5, to = 14, seq = "YIAKQRPISF",
                   D = c(1.1, 1.2, 1.15, 2.4, 2.5, 2.45, 0.9, 1.0, 0.95, 2.0, 2.1, 2.05))
d_long <- as_hdx_long(long, state = "state", time = "t_min", replicate = "run",
                      start = "from", end = "to", sequence = "seq", uptake = "D",
                      time_unit = "min")
output_tp(d_long)
#> Protein.States used: apo holo
#> Deut.times used: 60.00s 600.00s
#> Number of replicates used: apo=3 holo=3
#> no csv file written
#>   Deut.Time Start End   Sequence Search.RT Charge apo_X..Deut_1 apo_X..Deut_2
#> 1    60.00s     5  14 YIAKQRPISF         0      1           1.1           1.2
#> 2   600.00s     5  14 YIAKQRPISF         0      1           2.4           2.5
#>   apo_X..Deut_3 holo_X..Deut_1 holo_X..Deut_2 holo_X..Deut_3
#> 1          1.15            0.9            1.0           0.95
#> 2          2.45            2.0            2.1           2.05
```

## Replicates

All replicates are used, and states can have different numbers of
replicates. In `three_states.csv` the wild type has 4 replicates; the
mutants have 3 or 4, depending on the peptide. The table has one column
per replicate (4 per state here), and a missing replicate is `NA`.

``` r

a3 <- output_tp(d3)
#> Protein.States used: WT_pH7.5 D109H_pH7.5 R120G_ph7.5
#> Deut.times used: 4.00s 60.00s 1800.00s
#> 66 peptide row(s) dropped: 66 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: WT_pH7.5=4 D109H_pH7.5=4 R120G_ph7.5=4
#> no csv file written
attr(a3, "replicates")
#>    WT_pH7.5 D109H_pH7.5 R120G_ph7.5 
#>           4           4           4
## replicates with a value, per state, for the first peptides
sapply(c("WT", "D109H", "R120G"), function(s) rowSums(!is.na(a3[1:5, grep(s, colnames(a3))])))
#>   WT D109H R120G
#> 1  4     3     3
#> 2  4     3     3
#> 3  4     3     3
#> 4  4     3     3
#> 5  4     3     3
```

- A peptide is used when it has at least 2 replicates in every state
  (`min_replicates`). Dropped peptides are reported.
- `replicates = 3` uses at most 3 replicates per state.

## Point mutants

Peptides covering a mutation have a different sequence in the mutant. By
default only peptides with the same sequence are compared.
`seq_match = TRUE` matches peptides by position, so these peptides are
kept:

``` r

by_seq <- output_tp(d3, times = "60.00s")
#> Protein.States used: WT_pH7.5 D109H_pH7.5 R120G_ph7.5
#> 22 peptide row(s) dropped: 22 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: WT_pH7.5=4 D109H_pH7.5=4 R120G_ph7.5=3
#> no csv file written
by_pos <- output_tp(d3, times = "60.00s", seq_match = TRUE)
#> Protein.States used: WT_pH7.5 D109H_pH7.5 R120G_ph7.5
#> Number of replicates used: WT_pH7.5=4 D109H_pH7.5=4 R120G_ph7.5=3
#> no csv file written
c(by_sequence = nrow(by_seq), by_position = nrow(by_pos))
#> by_sequence by_position 
#>          40          50
```

[`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md)
and
[`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md)
do the same with `by = "position"`.

## Value scales

`percent` sets the values used everywhere:

| `percent` | Values |
|----|----|
| `FALSE` (default) | uptake in Da |
| `TRUE` | percent deuteration from the input (HDXExaminer `Deut %`) |
| `"per_residue"` | uptake per exchangeable residue (Da / [`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)) |
| `"FD"` | percent of the fully deuterated control |

Long peptides take up more deuterium in Da than short ones. `TRUE`,
`"per_residue"` and `"FD"` remove this, so peptide length does not
affect the result. `"FD"` also corrects for back-exchange. p-values are
the same on every scale.

``` r

a_fd <- output_tp(d3, times = "60.00s", percent = "FD")
#> Protein.States used: WT_pH7.5 D109H_pH7.5 R120G_ph7.5
#> 22 peptide row(s) dropped: 22 not measured in every state, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: WT_pH7.5=4 D109H_pH7.5=4 R120G_ph7.5=3
#> no csv file written
summary(unlist(a_fd[, 7:ncol(a_fd)]))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   27.12   57.00   74.59   70.70   86.24  105.68
```

With `"FD"`, peptides without an FD value stay in the table with `NA`.

## Peptide checks

[`peptide_qc()`](https://mkajano.github.io/HDXBoxeR/reference/peptide_qc.md)
flags missing controls, uptake of 0 at every time point, negative
uptake, and odd FD values. It does not remove anything.

``` r

qc <- peptide_qc(d3)
qc[qc$issues != "", c("Protein.State", "Start", "End", "issues")]
#> [1] Protein.State Start         End           issues       
#> <0 rows> (or 0-length row.names)
```

## How this connects to the rest of the package

| After loading | Use |
|----|----|
| Check the data | [`peptide_qc()`](https://mkajano.github.io/HDXBoxeR/reference/peptide_qc.md), [`general_info()`](https://mkajano.github.io/HDXBoxeR/reference/general_info.md), [`plot_peptide_coverage_map()`](https://mkajano.github.io/HDXBoxeR/reference/plot_peptide_coverage_map.md) ([`vignette("data_quality")`](https://mkajano.github.io/HDXBoxeR/articles/data_quality.md)) |
| Organize for analysis | [`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md) (time points), [`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md) (time courses) |
| Compare two states | [`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md) and the plots in [`vignette("HDXBoxeR")`](https://mkajano.github.io/HDXBoxeR/articles/HDXBoxeR.md) |
| Compare more states, regions, time courses | [`compare_states_tp()`](https://mkajano.github.io/HDXBoxeR/reference/compare_states_tp.md), [`compare_regions()`](https://mkajano.github.io/HDXBoxeR/reference/compare_regions.md), [`compare_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/compare_timecourse.md) ([`vignette("multistate")`](https://mkajano.github.io/HDXBoxeR/articles/multistate.md)) |
| Report | [`hdx_summary_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_summary_table.md), [`hdx_results_table()`](https://mkajano.github.io/HDXBoxeR/reference/hdx_results_table.md), [`render_report()`](https://mkajano.github.io/HDXBoxeR/reference/render_report.md) |
