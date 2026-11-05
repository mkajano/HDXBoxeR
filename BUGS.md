# HDXBoxeR — bug worklist (branch correlation_V2)

Scratch review notes. Not a git-tracked deliverable. No R files were
modified.

All findings below were verified by reading the code and, where marked
**\[reproduced\]**, by actually running the code against
`inst/extdata/All_results_table.csv` with R 4.2.1.

------------------------------------------------------------------------

## Severity 1 — crashes / silently wrong results

### \[R/boxer_input_to_long.R:41\] `all_data_long_input()` reads an undefined object `path`

**Bug:** the function argument is `filepath`, but the body reads `path`.
There is no `path` object anywhere in the package, so the function
throws on every call. **Evidence:**

``` r
all_data_long_input<-function(filepath){
  b1<-read.csv(path)                 # <- line 41: `path` is undefined
```

**\[reproduced\]**
`Error in read.table(file = file, ...) : object 'path' not found`
**Impact:** the function is 100% non-functional. Every call errors.
**Suggested fix:** `b1 <- read.csv(filepath)`.

------------------------------------------------------------------------

### \[R/boxer_input_to_long.R:22-24\] `percent = TRUE` mangles `Protein.State` and `Experiment`

**Bug:** the state/replicate split is done by splitting the column name
on the literal letter `"X"`. That only works for uptake columns
(`Unbound_X..Deut_1`). With `percent = TRUE`
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
produces `Unbound_Deut.._1`, which contains no `"X"`, so
`str_split_fixed(..., 'X', 2)` returns the whole string in column 1 and
`""` in column 2. **Evidence:**

``` r

nms<-data.frame(str_split_fixed(v5$Experiment, 'X', 2))
v5$Protein.State <- str_sub(nms[,1], end = -2)      # line 23
v5$Experiment<-paste("X", nms[,2], sep="")          # line 24
```

**\[reproduced\]** `boxer_input_to_long(file_nm, percent = TRUE)` gives
`Protein.State = "Unbound_Deut.._"` / `"bound_Deut.._"` and
`Experiment = "X"` for all 768 rows (replicate identity completely
lost). Correct output (`percent = FALSE`) is
`Protein.State = "Unbound"/"bound"`, `Experiment = "X..Deut_1..3"`.
**Impact:** `percent = TRUE` is a documented argument and silently
produces a long-format table with unusable state labels and no replicate
distinction. Any downstream group-by on `Protein.State`/`Experiment` is
wrong, no error raised. **Suggested fix:** don’t parse on `"X"`. Split
on the known suffix instead, e.g.
`sub("_(X\\.\\.Deut|Deut\\.\\.)_[0-9]+$", "", nm)` for the state and
[`regmatches()`](https://rdrr.io/r/base/regmatches.html)/[`sub()`](https://rdrr.io/r/base/grep.html)
for the replicate number, or better, carry the state/replicate mapping
out of
[`output_tp()`](https://mkajano.github.io/HDXBoxeR/reference/output_tp.md)
directly.

------------------------------------------------------------------------

### \[R/boxer_input_to_long.R:27 and :44\] `frac_deut` is nonsense when `percent = TRUE`

**Bug:** `frac_deut` divides the value column by the number of
exchangeable protons. That is only meaningful for uptake in Daltons.
With `percent = TRUE` the column holds percent deuteration, so the
result is “percent per residue”. **Evidence:**

``` r

v5$nb_prot<-nb_exch_deut(v5)
v5$frac_deut<-as.numeric(v5$X..Deut)/v5$nb_prot   # line 27
```

**\[reproduced\]** `percent = TRUE` gives
`frac_deut = 6.42, 5.91, 3.77, ...` (should be 0–1). `percent = FALSE`
correctly gives `0.63, 0.61, 0.32`. **Impact:** silently wrong numbers
in a column named “fraction deuterated”. **Suggested fix:** compute
`frac_deut` only when `percent = FALSE`, or divide by `100 * nb_prot` in
the percent branch — or simply
[`stop()`](https://rdrr.io/r/base/stop.html) for `percent = TRUE`.

------------------------------------------------------------------------

### \[R/output_prep.R:78\] `FD2` writes `100` into the **standard-deviation** column, and `100` as an uptake in Daltons

**Bug:** the diff changed three `rep(0, ...)` to `rep(100, ...)`. Only
one of those three columns survives into the output as a *mean*; another
survives as a *standard deviation*. The FD rows now report
`st.dev_Deut_Uptake = 100`. In addition, in the default
`percent = FALSE` branch the mean column is an uptake in **Daltons**, so
`100` is not “100 %” — it is 100 Da. **Evidence:**

``` r

FD2<-data.frame(FD[,c(1,2,4:6,8,7,9)], rep(100,times=dim(FD)[1]),
                rep(100,times=dim(FD)[1]), rep(100,times=dim(FD)[1]))   # line 78
...
colnames(FD2)<-nm_final   # cols 9,10,11 -> st.dev_Peptide_Mass, Mean.Deut.Uptake, st.dev_Deut_Uptake
out<-out[,c(1:8,10:11)]   # line 123 drops col 9, KEEPS 10 (mean) and 11 (sd)
```

**\[reproduced\]** FD rows in the written csv:

    Protein.State Deut.Time ... Mean.Deut.Uptake_.Da. st.dev_Deut_Uptake
          Unbound        FD ...                   100                100

and in `percent = TRUE` mode:
`Mean.Percent_Deut_.Da. = 100, st.dev_Percent_Deut = 100`. **Impact:**
[`output_prep()`](https://mkajano.github.io/HDXBoxeR/reference/output_prep.md)
is explicitly “output ready for publication” (Masson et al. 2019
format). Every FD row now ships a bogus SD of 100, and in uptake mode a
bogus mean of 100 Da. This is a *regression* introduced by the
uncommitted diff — before it the values were 0, which was wrong for
percent mode but at least not a fake SD. **Suggested fix:** make the
three constants branch-aware. Something like: SD column always `0`; mean
column `100` when `percent = TRUE` and `FD[,10]` (the real FD `X..Deut`)
when `percent = FALSE`. i.e. move the `FD2`/`undeut2` construction
inside the `if (percent==FALSE) ... else ...` block where `nm_final` is
already branch-specific.

------------------------------------------------------------------------

### \[R/sd_ci_timepoint.R:24-26\] CI→SD conversion is ~2.5× too large (missing factor of 2, wrong t quantile)

**Bug:** the cited source (Cochrane 7.7.3.2, in the roxygen block) gives
`SE = (UL - LL) / (2 * t_{1-alpha/2, df})` and `SD = SE * sqrt(N)`. The
code omits the factor of 2 **and** uses `qt(alpha, ...)` (a one-tailed
quantile) where the two-sided `t.test(conf.level = 1-alpha)` interval
requires `qt(alpha/2, ...)`. **Evidence:**

``` r

tvalue=abs(qt(alpha, replicates*2-2))                     # line 24: should be alpha/2
tt<-t.test(x1, x2, conf.level = c(1-alpha))               # two-sided interval
sd.ci<-(tt$conf.int[2]-tt$conf.int[1])/tvalue*sqrt(replicates)   # line 26: missing /2
```

For `alpha = 0.01`, `replicates = 3`: code uses
`t = |qt(0.01,4)| = 3.747`, correct is
`2 * |qt(0.005,4)| = 2 * 4.604 = 9.21`. Ratio = **2.46× overestimate**.
**Impact:**
[`sd_ci_timepoint()`](https://mkajano.github.io/HDXBoxeR/reference/sd_ci_timepoint.md)
is the **default** SD source for
[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
(`sd_source = "CI"`). Every uncertainty ellipse is ~2.5× too big and the
global CI band derived from those SDs (`CI_2pts`) is ~2.5× too wide, so
significant peptides are under-called. This is a statistical-correctness
bug, not cosmetic. **Suggested fix:**

``` r

tvalue <- abs(qt(alpha/2, replicates*2-2))
sd.ci  <- (tt$conf.int[2]-tt$conf.int[1]) / (2*tvalue) * sqrt(replicates)
```

(Also consider using `tt$parameter`, the Welch df, rather than
`replicates*2-2`.) Update the formula in the roxygen block and in
`plot_fissure_vent.R`’s `@details`, which repeat the same wrong formula.

------------------------------------------------------------------------

### \[R/plot_fissure_vent.R:123\] The plot’s y axis is built from **state1**, not state2

**Bug:** both arguments to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) are
`av_state1`. The y-axis range is therefore derived from state1’s data,
not state2’s. **Evidence:**

``` r

plot(av1[, av_state1], av1[, av_state1], type = "n", xlab = "", ylab = "",
     xlim = xlim, ylim = ylim, ...)                       # line 123
```

**Impact:** whenever the two states have different value ranges (which
is the entire point of the plot) ellipses for state2 are drawn partly or
wholly off-canvas, with no warning. Inherited verbatim from
`plot_cor_pts_elps` (`R/correlation_functions_2states.R:24` has the
identical `plot(av1[,7], av1[,7], ...)`), so it was *not* fixed by
868a44f. **Suggested fix:**
`plot(av1[, av_state1], av1[, av_state2], ...)`.

------------------------------------------------------------------------

### \[R/residue_deconvolution.R:76-82\] `peptide_exch_residues()` includes one residue too many — contradicts `nb_exch_deut()`

**Bug:** the roxygen says the helper matches
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)’s
convention (“the first residue of each peptide and all Prolines are
excluded”).
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
actually excludes the first **two** residues:
`(End-Start) - nP + firstP - 1` = `L - 2 - nP_after_first`. The helper
starts at local index 2, giving `L - 1 - nP_after_first` — one extra
residue per peptide. **Evidence:**

``` r

peptide_exch_residues <- function(start, seq_string) {
  chars <- strsplit(seq_string, "")[[1]]
  if (length(chars) < 2) return(integer(0))
  local_idx <- 2:length(chars)                 # <- line 79, should be 3:length(chars)
  resid_idx <- start + local_idx - 1
  resid_idx[chars[local_idx] != "P"]
}
```

vs `R/nb_exchangeable_protons.R:15-16`:

``` r

nb_exch_deut<-c(c(df$End-df$Start)-str_count(df$Sequence,'P') ) +
  as.numeric(substr(df$Sequence, 1, 1)=="P")-1
```

**\[reproduced\]** for all 128 peptides in the shipped example,
`length(peptide_exch_residues(...)) - nb_exch_deut(...) == 1` (every
single row). **Impact:** the design matrix `A` spreads each peptide’s
observed uptake over one extra residue, so every per-residue estimate is
systematically diluted at the peptide’s N-terminal shoulder, and residue
`Start+1` (which the rest of the package treats as non-exchangeable) is
assigned a spurious value. Also breaks the vignette’s cross-validation
section, which copy-pastes the same helper (`vignettes/HDXBoxeR.Rmd`,
“Choosing `lambda_smooth`”) — so the CV is fit and scored with the same
off-by-one and cannot detect it. **Suggested fix:**

``` r

if (length(chars) < 3) return(integer(0))
local_idx <- 3:length(chars)
```

and mirror the change in the vignette chunk. Better: factor the helper
out as a real (possibly non-exported) function so
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md),
the deconvolution, and the vignette all share one definition —
`length(peptide_exch_residues(...)) == nb_exch_deut(...)` would then be
a natural unit test.

------------------------------------------------------------------------

## Severity 2 — inconsistent / misleading results, R CMD check breakage

### \[R/residue_level_average.R:58\] Coverage window is off by one and ignores Prolines

**Bug:** the comment claims the window matches “the exchangeable-residue
convention used elsewhere in this package”, but it uses `Start+1 .. End`
(excludes only residue 1) and does no Proline exclusion at all.
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
excludes the first two residues and all non-first Prolines;
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)
(after fixing the bug above) excludes first-two + Prolines.
**Evidence:**

``` r

## a residue is covered if it falls within (Start+1)..End, matching
## the exchangeable-residue convention used elsewhere in this package
covering <- (sub$Start + 1) <= r & sub$End >= r        # line 58
```

**\[reproduced\]** with the example data, residue 41 (a Proline — first
peptide is `GPLGSK...` at Start 40) is reported as covered by
[`residue_level_average()`](https://mkajano.github.io/HDXBoxeR/reference/residue_level_average.md)
but is `NA` in
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md).
**Impact:** the two functions are documented as directly comparable
baselines, but they disagree about which residues even exist. Prolines
get fabricated uptake values in the averaged output. **Suggested fix:**
use the same shared `peptide_exch_residues()` helper as
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)
to decide coverage, rather than an inline range test.

------------------------------------------------------------------------

### \[R/boxer_input_to_long.R:20,22\] `gather()` and `str_split_fixed()` are used but not imported

**Bug:** neither
[`tidyr::gather`](https://tidyr.tidyverse.org/reference/gather.html) nor
[`stringr::str_split_fixed`](https://stringr.tidyverse.org/reference/str_split.html)
appears in NAMESPACE. NAMESPACE only has `importFrom(tidyr, separate)`
and `importFrom(stringr, str_count)` / `importFrom(stringr, str_sub)`.
**Evidence:** `grep '^import' NAMESPACE` → no `gather`, no
`str_split_fixed`. **Impact:** `R CMD check` NOTE “no visible global
function definition for ‘gather’ / ‘str_split_fixed’”; at runtime the
function errors in an installed package unless the user happens to have
`tidyr`/`stringr` attached. It only worked in my test because I
[`library()`](https://rdrr.io/r/base/library.html)d them manually.
**Suggested fix:** add `@importFrom tidyr gather` and
`@importFrom stringr str_split_fixed` (and re-roxygenize). Note
[`tidyr::gather()`](https://tidyr.tidyverse.org/reference/gather.html)
is retired-superseded;
[`tidyr::pivot_longer()`](https://tidyr.tidyverse.org/reference/pivot_longer.html)
would be the forward-looking choice.

------------------------------------------------------------------------

### \[NAMESPACE\] `boxer_input_to_long` and `all_data_long_input` have man pages but are not exported

**Bug:** the NAMESPACE diff adds only `plot_fissure_vent`,
`residue_deconvolution`, `residue_level_average`. Neither function in
`boxer_input_to_long.R` carries an `@export` tag, yet
`man/boxer_input_to_long.Rd` and `man/all_data_long_input.Rd` were
generated. **Evidence:** `grep boxer_input_to_long NAMESPACE` → no
match; `man/boxer_input_to_long.Rd` exists with a full `\usage{}` block.
**Impact:** users cannot call
[`HDXBoxeR::boxer_input_to_long()`](https://mkajano.github.io/HDXBoxeR/reference/boxer_input_to_long.md).
`R CMD check` also emits “Rd sections for objects not exported”. (Four
pre-existing Rd files have the same problem:
`deuteration_woods_timecourse_frac`, `is.nan.data.frame`,
`plot_cor_pts_elps`, `sd_ci_timepoint` — worth cleaning up in the same
pass.) **Suggested fix:** add `@export` to both roxygen blocks (they
clearly document user-facing arguments), or add `@keywords internal` /
`@noRd` if they are meant to be internal.

------------------------------------------------------------------------

### \[tests/ — missing `tests/testthat.R`\] the test suite never runs

**Bug:** `tests/testthat/test-residue_deconvolution.R` exists but there
is no `tests/testthat.R` driver file. DESCRIPTION correctly has
`testthat (>= 3.0.0)` in Suggests and `Config/testthat/edition: 3`, but
without the driver `R CMD check` finds no runnable test file in
`tests/`. **Evidence:** `ls tests tests/testthat` → `tests/` contains
only the `testthat/` directory. **Impact:** all four tests are dead
weight — `R CMD check` reports 0 tests run, and CI would pass on a
completely broken
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md).
**Suggested fix:** create `tests/testthat.R` containing:

``` r

library(testthat)
library(HDXBoxeR)

test_check("HDXBoxeR")
```

(`usethis::use_testthat(3)` does this.)

**Separately, on the tests themselves:** they run and pass conceptually
(I checked the assertions against real output: `nrow == 225*4 == 900` ✓;
the two `expect_error` cases both really do hit the
[`stop()`](https://rdrr.io/r/base/stop.html) because the two column-name
regexes are mutually exclusive — `_X\.\.Deut_[0-9]+$` vs
`_Deut\.\._[0-9]+$` ✓; deconvolution values fall in `[0,1]` ✓; averaged
values fall in `[10.8, 102.8]` ✓). But every assertion is a shape/bounds
check. Nothing tests numerical correctness, which is why the off-by-one
in `peptide_exch_residues()` (bug above) sails straight through. Worth
adding: - a synthetic-data round trip: build peptides with known
per-residue uptake, check `A %*% x_hat ≈ b`; -
`length(peptide_exch_residues(Start, Sequence)) == nb_exch_deut(row)`
for every row; - `lambda_smooth = 0` vs `> 0` actually changes the
result (currently untested); - an `NA` residue for a region with no
peptide coverage.

------------------------------------------------------------------------

### \[man/plot_fissure_vent.Rd\] Rd file is stale relative to the R source — `R CMD check` WARNING

**Bug:** the `.Rd` was generated from an older version of the roxygen
block. It documents an argument `sd` that does not exist, and omits four
that do. **Evidence:** `\usage{}` lists `position_col`, `deut`,
`length1`, `xlim`; `\arguments{}` has no `\item` for any of them, but
does have `\item{sd}{flag. CI (default) ...}` — the argument is actually
named `sd_source`. The `\examples{}` block calls
[`woods_CI_plot_frac()`](https://mkajano.github.io/HDXBoxeR/reference/woods_CI_plot_frac.md)
and
[`output_tc()`](https://mkajano.github.io/HDXBoxeR/reference/output_tc.md),
i.e. it never calls
[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
at all. **Impact:** `R CMD check` → WARNING “Undocumented arguments in
documentation object ‘plot_fissure_vent’: ‘position_col’ ‘deut’
‘length1’ ‘xlim’” and “Documented arguments not in … : ‘sd’”. The R
source is fine (it has `@param sd_source`, `@param xlim`); only the Rd
is stale. But `position_col`, `deut` and `length1` are undocumented in
the source too. **Suggested fix:** re-run `roxygen2::roxygenise()`, and
add `@param position_col`, `@param length1` to `R/plot_fissure_vent.R`.
`deut` (line 55) is accepted but never referenced anywhere in the body —
delete it.

------------------------------------------------------------------------

### \[R/plot_fissure_vent.R:60\] Default `state1`/`state2` are truncated by an off-by-one `str_sub()`

**Bug:** `end = -10` strips 9 trailing characters, but the suffix on an
averaged column (`av_Unbound_X..Deut`) is `_X..Deut`, 8 characters.
**Evidence:**

``` r

nm1 <- str_sub(colnames(av1)[7:ncol(av1)], start = 4, end = -10)   # line 60
if (missing(state1)) state1 <- nm1[1]
```

**\[reproduced\]** uptake input → `nm1 = c("Unboun", "boun")`; percent
input (suffix `_Deut..`, 7 chars) → `nm1 = c("Unbou", "bou")`.
**Impact:** currently masked, because
`resolve_col()`/`resolve_raw_cols()` use
[`startsWith()`](https://rdrr.io/r/base/startsWith.html) and a truncated
prefix still matches. It stops being masked as soon as two states share
a prefix that only differs in the last character(s) — e.g. `"State1"` /
`"State2"` truncate to `"State"` / `"Stat"`, and `resolve_col()` then
throws “State ‘State’ matches multiple columns”. Also the truncated
string is what appears in any user-facing error message. **Suggested
fix:** don’t count characters. Derive the root by stripping the known
suffix:
`sub("^av_", "", sub("_(X\\.\\.Deut|Deut\\.\\.)$", "", colnames(av1)[7:ncol(av1)]))`.

------------------------------------------------------------------------

### \[R/plot_fissure_vent.R:49\] The `@examples` block passes a numeric `times` and silently produces an empty plot

**Bug:** `Deut.Time` values are character (`"3.00s"`, `"60.00s"`, …).
The example passes `times = 3.00`.
[`select_indices()`](https://mkajano.github.io/HDXBoxeR/reference/select_indices.md)
compares `df$Deut.Time == 3.00`, which is `FALSE` everywhere.
**Evidence:**

``` r

#' plot_fissure_vent(df = a, replicates = 3, times = 3.00, alpha = 0.01,
#'                    pvalue = 0.01, state1 = "Unbound", state2 = "bound")   # line 49-50
```

**\[reproduced\]** `times = 3.00` → `ind1` empty →
[`CI_2pts()`](https://mkajano.github.io/HDXBoxeR/reference/CI_2pts.md)
on a zero-length vector → `message: "time 3: CI = NaN"`, plot drawn with
zero peptides, no error or warning. `times = "3.00s"` →
`"time 3.00s: CI = 1.607"` and peptides drawn. **Impact:** the
documented example produces a blank plot. More importantly the function
gives no feedback for an unmatched `times` value, so a user typo yields
a silently empty figure. **Suggested fix:** change the example to
`times = "3.00s"`, and add a guard in the loop:
`if (length(ind1) == 0) { warning(sprintf("no peptides match times = '%s'", times[i])); next }`.

------------------------------------------------------------------------

### \[R/robot_indexes_df.R:43-44\] `replicates` argument is ignored when computing the CI

**Bug:** `replicates = 3` is hardcoded in the
[`prep_timecourse_plot_sd()`](https://mkajano.github.io/HDXBoxeR/reference/prep_timecourse_plot_sd.md)
call, while the sibling function
[`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md)
(line 99-100 of the same file) correctly passes
`replicates = replicates`. **Evidence:**

``` r

# robot_indexes_df, line 43:
CI_all<-prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates=3, alpha=alpha)
# robot_indexes, line 99:
CI_all<-prep_timecourse_plot_sd(control_df_up, variant_df_up, replicates=replicates, alpha)
```

**Impact:** for any dataset with ≠3 replicates,
[`robot_indexes_df()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes_df.md)
computes the CI threshold with the wrong t-value and wrong pooled-SD
normalization, so it returns the wrong significant-peptide set — while
[`pv_timecourse()`](https://mkajano.github.io/HDXBoxeR/reference/pv_timecourse.md)
and
[`prep_timecourse_plot_ave()`](https://mkajano.github.io/HDXBoxeR/reference/prep_timecourse_plot_ave.md)
in the same function use the correct `replicates`. Silent. **Suggested
fix:** `replicates = replicates`.

------------------------------------------------------------------------

### \[R/robot_indexes_df.R:26-56 and :82-109\] Multi-state loop discards all but the last state

**Bug:** `peptide_all <- c()` is re-initialized inside the
`for (state in ...)` loop, and `sig_pep` is built *after* the loop. With
more than 2 protein states, the significant peptides found for every
state except the last are thrown away. **Evidence:**

``` r

for ( state in states[2:length(states)]) {
    ...
    peptide_all<-c()                                   # line 46: reset each iteration
    for (i in 7:dim(sh_avc)[2]) { peptide_all<-c(peptide_all, which(...)) }
    peptide_all<-sort(unique(peptide_all))
  }                                                    # line 52: loop ends
  sig_pep<-data.frame(1:length(peptide_all), sh_avc[peptide_all,c(1:3,5)])   # line 54
```

Identical structure in
[`robot_indexes()`](https://mkajano.github.io/HDXBoxeR/reference/robot_indexes.md)
(lines 82-109). **Impact:** the docs say “As default all states are
chosen”, but the returned peptide set reflects only `states[1]` vs
`states[length(states)]`. Also note `sh_avc` used on line 54 is the last
iteration’s version — consistent by accident, not by design. **Suggested
fix:** move `peptide_all <- c()` above the `for (state ...)` loop and
drop the mid-loop reset (or accumulate into a list keyed by state and
return that).

------------------------------------------------------------------------

### \[R/correlation_functions_2states.R:38-39\] The ellipse bug fixed in `plot_fissure_vent` still exists here

**Bug:**
[`plot_cor_pts_elps()`](https://mkajano.github.io/HDXBoxeR/reference/plot_cor_pts_elps.md)
is the function
[`plot_fissure_vent()`](https://mkajano.github.io/HDXBoxeR/reference/plot_fissure_vent.md)
was derived from. Its first branch divides the ellipse **y-centre** by 2
and applies `/2` to only one of the two radii — exactly the bug commit
868a44f fixed in `plot_fissure_vent.R`. The second branch (line 51) has
it right. **Evidence:**

``` r

# line 38-39, is.na(position_col) branch  -- WRONG
DrawEllipse(av1[i,state1], av1[i,state2]/2, sd1[i,state1]/2,
            sd1[i,state2], ...)
# line 51-52, else branch                 -- CORRECT
DrawEllipse(av1[i,state1], av1[i,state2], sd1[i,state1]/2, sd1[i,state2]/2,
            col = position_col[i], border= position_col[i], ...)
```

**Impact:** the default call path of
[`plot_cor_pts_elps()`](https://mkajano.github.io/HDXBoxeR/reference/plot_cor_pts_elps.md)
(no `position_col` supplied) draws every ellipse at half its correct y
position with a y-radius 2× too large. `plot_cor_pts_elps` is used by
the exported `plot_correlation_functions` family. **Suggested fix:**
make line 38 match line 51:
`DrawEllipse(av1[i,state1], av1[i,state2], sd1[i,state1]/2, sd1[i,state2]/2, ...)`.
Same file line 24 also has the `plot(av1[,7], av1[,7], ...)`
duplicated-axis bug.

------------------------------------------------------------------------

## Severity 3 — edge cases, robustness

### \[R/plot_fissure_vent.R:140,155-156\] Significance threshold is `CI/2`, but the rest of the package uses `CI`

**Bug:** the roxygen block claims the peptide outline uses “the same
two-part standard as the rest of this package (CI_2pts + a p-value
cutoff)”. It doesn’t — it halves the CI. **Evidence:**

``` r

is_sig <- (delta > ci / 2) && !is.na(pv_state1v2[pep]) && (pv_state1v2[pep] < pvalue)  # line 140
points(c(-10, 110), c(-10, 110) + ci / 2, type = "l", col = "grey")                    # line 155
```

vs everywhere else in the package:

``` r
R/inside_functions_plot.R:45   abs1 <- abs(df_av[,i]) > CI[i-7]
R/list_significantly_different_peptides.R:33  abs.a <- abs(df1[...,i]) > CI
R/robot_indexes_df.R:49        abs(sh_avc_up[,i]-sh_avv_up[,i]) > CI_all[i-6]
```

**Impact:** twice as many peptides are flagged significant here as by
`list_significantly_different_peptides()` / the robot plots on the same
data. The `/2` was inherited from `plot_cor_pts_elps` (where it was only
a drawing convention with no significance call attached);
`plot_fissure_vent` promoted it to a statistical decision. Note the
~2.5× SD inflation in `sd_ci_timepoint` (Severity 1) partly cancels
this, which is why it doesn’t look obviously wrong on the example data.
**Suggested fix:** decide which convention is intended and make band +
test agree — either draw at `±ci` and test `delta > ci`, or keep the
visual `±ci/2` band but document clearly that it is a half-band and test
against `ci`.

------------------------------------------------------------------------

### \[R/residue_level_average.R:51\] `1 / nb_exch_deut()` can be `Inf` or negative

**Bug:** no guard on the weight denominator.
[`nb_exch_deut()`](https://mkajano.github.io/HDXBoxeR/reference/nb_exch_deut.md)
returns `L - 2 - nP_after_first`, which is `0` for a 2-residue peptide
(or e.g. `PGP`-heavy short peptides) and can go negative for a short
Proline-rich peptide. **Evidence:**

``` r

w <- 1 / nb_exch_deut(sub)     # line 51
...
res_mat[r, ] <- colSums(vals * wi * ok, na.rm = TRUE) / colSums(wi * ok, na.rm = TRUE)
```

**Impact:** a zero → `Inf` weight makes the whole residue
`Inf/Inf = NaN`; a negative weight silently produces a weighted mean
outside the data range. Not triggered by the shipped example
(`range(nb_exch_deut(a)) == c(5, 19)`), so it will show up only on a
user’s dataset with short peptides — the worst kind of latent bug.
**Suggested fix:** drop rows with `nb_exch_deut(sub) <= 0` before
weighting, and warn.

------------------------------------------------------------------------

### \[R/residue_level_average.R:63\] Uncovered/all-NA columns yield `NaN`, not `NA`

**Bug:** if every covering peptide has `NA` in a given column, the
denominator `colSums(wi * ok, na.rm = TRUE)` is `0` and the result is
`NaN` rather than the documented `NA`. **Evidence:**

``` r

res_mat[r, ] <- colSums(vals * wi * ok, na.rm = TRUE) / colSums(wi * ok, na.rm = TRUE)  # line 63
```

Doc says: “`NA` where no peptide covers a residue.” **Impact:**
[`is.na()`](https://rdrr.io/r/base/NA.html) catches `NaN` so most
downstream code is fine, but plotting and
[`identical()`](https://rdrr.io/r/base/identical.html)/comparison code
that distinguishes them will not be. Minor. **Suggested fix:**
`den <- colSums(wi*ok, na.rm=TRUE); res_mat[r, ] <- ifelse(den > 0, num/den, NA_real_)`.

------------------------------------------------------------------------

### \[R/sd_ave_diff_pv_tp.R:17,47,82 and R/sd_ci_timepoint.R:17\] `str_sub(nm_root, end = -3)` breaks with ≥10 replicates

**Bug:** the root name is derived by chopping a fixed 2 characters,
assuming the replicate suffix is always `_1`..`_9`. **Evidence:**
(identical line in all four functions)

``` r

nm_root<-colnames(df[,c(6+((1:nb_sets)-1)*replicates+1)])
nm_root<-str_sub(nm_root, end = -3)      # "Unbound_X..Deut_1" -> "Unbound_X..Deut"
```

**Impact:** with `replicates >= 10` the first column of a set is still
`_1`, so this happens to survive — but any code path that reaches `_10`
gets `"Unbound_X..Deut_"`. Combined with the `end = -10` truncation in
`plot_fissure_vent.R:60`, column-name parsing in this package is
character-count based throughout and will misbehave on unusual state
names. Low priority but worth a systematic fix. **Suggested fix:**
`sub("_[0-9]+$", "", nm_root)`.

------------------------------------------------------------------------

### \[R/plot_fissure_vent.R:56,123,145\] `...` is forwarded to both `plot()` and `DrawEllipse()`

**Bug:** `@param ...` is documented as “other plotting parameters,
passed to DrawEllipse”, but the same `...` also goes to
[`plot()`](https://rdrr.io/r/graphics/plot.default.html). **Evidence:**

``` r

plot(av1[, av_state1], av1[, av_state1], type = "n", ..., ...)   # line 123
DrawEllipse(..., col = fill_col, border = border_col, ...)       # line 145
```

**Impact:** any `DrawEllipse`-specific argument (`nv`, `plot`, `angle`)
reaches [`plot.default()`](https://rdrr.io/r/graphics/plot.default.html)
→ [`par()`](https://rdrr.io/r/graphics/par.html) and produces “… is not
a graphical parameter” warnings on every call; conversely `plot`-only
arguments reach `DrawEllipse`. Also `col`/`border` are set explicitly on
line 145, so a user passing `col =` gets “formal argument matched by
multiple actual arguments”. **Suggested fix:** either restrict `...` to
`DrawEllipse` (and give
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) its own
explicit args) or add a separate `ellipse_args = list()` parameter.

------------------------------------------------------------------------

### \[vignettes/HDXBoxeR.Rmd, “Choosing `lambda_smooth`”\] The CV chunk is `eval`’d and re-fits the deconvolution 12 times

**Bug:** the cross-validation chunk has no `eval=FALSE` / `cache=TRUE`,
and runs
[`residue_deconvolution()`](https://mkajano.github.io/HDXBoxeR/reference/residue_deconvolution.md)
for 4 λ values × 3 folds. **Evidence:**
`for (lam in lambda_grid) { for (f in 1:k) { fit <- residue_deconvolution(train_df, protein_length = 225, lambda_smooth = lam) ... } }`
**Impact:** adds ~12× a single fit (≈0.6 s each on the shipped
128-peptide example, so tolerable here) to every `R CMD check` /
vignette build, and scales badly for anyone who adapts it. It also
silently re-defines `peptide_exch_residues()` with the same off-by-one
as the package internal, so the CV cannot detect that bug. **Suggested
fix:** `eval = FALSE` on the chunk (it’s illustrative), and call the
shared package helper instead of re-defining it.
