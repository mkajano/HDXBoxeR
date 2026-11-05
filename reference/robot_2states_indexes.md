# Returns a robot plot for selected peptides for 2 protein states.

Modification of butterfly plot. x axis residues. y axis % deuteration
for one variant above the axis and for second peptide below the axis.
Peptides are compared between the sets for the significance change
between sets. If there is significant change beteween sets peptides are
plotted for all timepoints. Significanty different timepoints for the
peptides are colored. Peptides ranges are plotted as a line at
corresponding % deuteration values.

## Usage

``` r
robot_2states_indexes(
  thP,
  th,
  indexes,
  states,
  replicates = 3,
  pvalue = 0.01,
  alpha = 0.01,
  ylim,
  xlim,
  p_adjust = "none"
)
```

## Arguments

- thP:

  output of output_tc(percent = TRUE). Raw data for procent deuteration
  for time courses

- th:

  output of output_tc(). Raw data for uptake deuteration for time
  courses

- indexes:

  indexes of peptides to be drawn.

- states:

  Need to choose only two protein states

- replicates:

  number of replicates in sample. Default set to 3.

- pvalue:

  p-value cutoff here set up to 0.01

- alpha:

  critical interval, to have more restrictive use lower values,
  default=0.01

- ylim:

  y-axis range

- xlim:

  x-axis range. Set as default from max and minimum residues for the
  protein

- p_adjust:

  multiple-testing correction of the Welch p-values, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html): `"none"`
  (default, raw p-values) or e.g. `"BH"`. A peptide is significant when
  its (adjusted) p-value is below the cutoff AND its difference exceeds
  the critical interval.

## Value

Robot maps for timecourses for 2 protein states and selected indexes.

## Examples

``` r
file_nm<-system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
tm_df<-output_tc(filepath=file_nm)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
tmP_df<-output_tc(filepath=file_nm, percent=TRUE)
#> Protein.States used: Unbound bound
#> Deut.times used: 3.00s 60.00s 1800.00s 72000.00s
#> 1 peptide row(s) dropped: 1 not measured in every time point, 0 with fewer than 2 replicate value(s)
#> Number of replicates used: 3.00s=3 60.00s=3 1800.00s=3 72000.00s=3
#> no csv file written
names_states<- nm_states(file_nm) ### returns states names
ind1<-robot_indexes(thP = tmP_df, th=tm_df, pvalue=0.001, alpha=0.01, states=names_states[1:2])
robot_2states_indexes(thP = tmP_df, th=tm_df,
 states=names_states[1:2],indexes =ind1, pvalue=0.001, alpha=0.001)
```
