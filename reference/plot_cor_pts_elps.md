# Correlation functions plots and prep data

Correlation of two states, the st.dev are shown as ellipsoids.

## Usage

``` r
plot_cor_pts_elps(
  df,
  replicates = 3,
  times,
  state1,
  state2,
  position_col = NA,
  length1 = NA,
  xlim = NULL,
  ylim = NULL,
  ...
)
```

## Arguments

- df:

  input from tp function

- replicates:

  number of replicates in sample. Default set to 3.

- times:

  times to be used

- state1:

  First state

- state2:

  second state

- position_col:

  default NA, flag that allows to color based on domain, if TRUE then
  the length1 vector needs to be provided

- length1:

  vector of numeric

- xlim:

  x-axis range

- ylim:

  y-axis range

- ...:

  other plotting parameters

## Value

Correlation plot
