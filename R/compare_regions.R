#' Region-level ANOVA: compare protein states within residue regions
#'
#' Tests whether protein states differ in deuterium uptake within each
#' region you define. Each peptide acts as its own control: the state
#' effect is a two-way ANOVA `uptake ~ peptide + state` on the per-peptide
#' mean uptake (replicates averaged), so it is the uptake difference
#' *within* peptides, averaged over the region's peptides, judged against
#' how consistent that difference is from peptide to peptide (peptides are
#' the replication unit, as in a randomized-block / repeated-measures
#' design; with two states this is a paired t-test of the per-peptide
#' differences). This is far more powerful than pooling peptide means into
#' one distribution per state (as [compare_states_unmatched()] must do when
#' peptide maps differ), because differences between peptides (length,
#' location) are removed before the states are compared.
#'
#' Requires the same peptides in every compared state (matched peptides,
#' as for [output_tp()]); peptides missing from any state are left out of
#' that region. For states with different peptide maps use
#' [compare_states_unmatched()].
#'
#' For each region (and time point) the result reports:
#' \itemize{
#'   \item an F-test of the state effect (do the states differ at all in
#'     this region; with two states it is the test of their difference),
#'   \item an F-test of the state x peptide interaction on the replicate
#'     values (does the state effect differ between peptides of the region
#'     by more than replicate noise - a significant value means the effect
#'     is localized or mixed in direction, so the average shift alone does
#'     not describe the region; look at the peptides with
#'     [compare_states_tp()] or [plot_fissure_vent()]),
#'   \item pairwise state contrasts: average shift (Da, state2 - state1) per
#'     peptide with its standard error and t-test p-value from the model.
#' }
#' With several time points and `time_model = TRUE`, one model per region
#' is fitted over all of them on the per-peptide-and-time means,
#' `uptake ~ peptide * time + state`, testing the state effect averaged over
#' time, plus whether it changes with time (state x time).
#'
#' A region needs at least two peptides. Every peptide weighs equally,
#' whatever its replicate noise. p-values are raw by default; `p_adjust`
#' corrects across all regions (and time points), and separately across
#' all pairwise contrasts.
#'
#' @param filepath input file or table, as for [output_tp()].
#' @param times one or more `Deut.Time` values. Default: all.
#' @param regions named list of `c(start, end)` residue ranges, e.g.
#'   `list(NTR = c(1, 85), ACD = c(86, 165))`. A peptide belongs to a region
#'   when it lies entirely inside it.
#' @param states states to compare (two or more). Default: all.
#' @param reference optional state; pairwise contrasts are then only
#'   `reference` vs each other state.
#' @param by peptide matching between states: `"sequence"` (default) or
#'   `"position"`, as in [compare_states_tp()].
#' @param replicates optional maximum number of replicates per state.
#' @param percent value scale, as in [output_tp()]: `FALSE` (Da, default),
#'   `TRUE` or `"FD"`.
#' @param time_model `FALSE` (default): one model per region and time point.
#'   `TRUE`: one model per region over all `times` (needs two or more).
#' @param p_adjust multiple-testing correction passed to [stats::p.adjust()]:
#'   `"none"` (default) or e.g. `"BH"`.
#' @return a list of two data frames. `summary`: one row per region (and
#'   time point): `region`, `Deut.Time` (`"all"` with `time_model = TRUE`),
#'   `n_peptides`, `n_obs`, `F_state`, `df_state`, `df_resid`, `p_state`,
#'   `p_state_adj`, `p_state_x_peptide`, and with `time_model = TRUE` also
#'   `p_state_x_time`. `pairs`: one row per region, time and state pair:
#'   `estimate` (Da, or percent, state2 - state1), `se`, `p_value`, `p_adj`.
#' @seealso [compare_states_tp()] (peptide by peptide),
#'   [compare_states_unmatched()] (different peptide maps),
#'   [residue_regions()] (regions found from the data)
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' res <- compare_regions(file_nm, times = "1800.00s",
#'                        regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' res$summary
#' res$pairs
#' @export
compare_regions <- function(filepath, times, regions, states = NULL, reference = NULL,
                            by = c("sequence", "position"), replicates = NULL,
                            percent = FALSE, time_model = FALSE, p_adjust = "none") {
  by <- match.arg(by)
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  if (missing(regions) || length(regions) == 0 || is.null(names(regions))) {
    stop("regions is required: a named list of c(start, end) residue ranges", call. = FALSE)
  }
  raw <- .hdx_raw(filepath)
  if (is.null(states)) states <- unique(arg_df(raw)$Protein.State)
  if (length(states) < 2) stop("need at least two states", call. = FALSE)
  if (!is.null(reference) && !reference %in% states) {
    stop("reference '", reference, "' not among states", call. = FALSE)
  }
  a <- suppressMessages(output_tp(raw, states = states, replicates = replicates,
                                  seq_match = (by == "position"), percent = percent))
  if (missing(times)) times <- unique(a$Deut.Time)
  a <- a[a$Deut.Time %in% times, , drop = FALSE]
  if (time_model && length(unique(a$Deut.Time)) < 2) {
    stop("time_model = TRUE needs two or more time points", call. = FALSE)
  }

  ## long format: one row per replicate value
  long <- do.call(rbind, lapply(states, function(s) {
    cols <- .hdx_state_cols(s, colnames(a))
    data.frame(Deut.Time = rep(a$Deut.Time, length(cols)),
               Start = rep(a$Start, length(cols)), End = rep(a$End, length(cols)),
               peptide = rep(paste(a$Start, a$End, a$Charge, a$Sequence, sep = "_"), length(cols)),
               state = s, value = as.numeric(unlist(a[, cols], use.names = FALSE)),
               stringsAsFactors = FALSE)
  }))
  long <- long[!is.na(long$value), , drop = FALSE]
  long$state <- factor(long$state, levels = states)

  pairs <- if (!is.null(reference)) {
    lapply(setdiff(states, reference), function(s) c(reference, s))
  } else utils::combn(states, 2, simplify = FALSE)

  fit_one <- function(d, rname, tlabel, with_time) {
    d$peptide <- factor(d$peptide)
    n_pep <- nlevels(d$peptide)
    empty <- list(summary = NULL, pairs = NULL)
    if (n_pep < 2 || length(unique(d$state)) < 2) {
      if (n_pep == 1) message("region '", rname, "' (", tlabel, "): only one peptide - not testable as a region")
      return(empty)
    }
    ## State effect: tested on per-peptide means, with peptides as the
    ## replication unit (randomized block / repeated measures). Testing it
    ## against replicate noise instead would overstate significance whenever
    ## the effect differs between peptides, which in HDX it nearly always does.
    if (with_time) d$time <- factor(d$Deut.Time)
    grp <- if (with_time) c("peptide", "time", "state") else c("peptide", "state")
    m <- stats::aggregate(d["value"], d[grp], mean)
    if (with_time) {
      main <- stats::lm(value ~ peptide * time + state, data = m)
      null <- stats::lm(value ~ peptide * time, data = m)
      p_st_time <- .anova_p(main, stats::lm(value ~ peptide * time + state + state:time, data = m))
      ## heterogeneity between peptides, against replicate noise
      p_het <- .anova_p(stats::lm(value ~ peptide * time + state * time, data = d),
                        stats::lm(value ~ peptide * time + state * time + state:peptide, data = d))
    } else {
      main <- stats::lm(value ~ peptide + state, data = m)
      null <- stats::lm(value ~ peptide, data = m)
      p_het <- .anova_p(stats::lm(value ~ peptide + state, data = d),
                        stats::lm(value ~ peptide * state, data = d))
    }
    av <- stats::anova(null, main)
    summ <- data.frame(region = rname, Deut.Time = tlabel, n_peptides = n_pep, n_obs = nrow(d),
                       F_state = av$F[2], df_state = av$Df[2], df_resid = av$Res.Df[2],
                       p_state = av$`Pr(>F)`[2], p_state_x_peptide = p_het,
                       stringsAsFactors = FALSE)
    if (with_time) summ$p_state_x_time <- p_st_time

    ## pairwise contrasts on the state coefficients (treatment coding)
    b <- stats::coef(main); V <- stats::vcov(main)
    lev <- levels(droplevels(m$state))
    coef_of <- function(s) if (s == lev[1]) 0 else b[[paste0("state", s)]]
    pr <- do.call(rbind, lapply(pairs, function(p) {
      if (!all(p %in% lev)) return(NULL)
      cvec <- stats::setNames(numeric(length(b)), names(b))
      if (p[2] != lev[1]) cvec[paste0("state", p[2])] <- 1
      if (p[1] != lev[1]) cvec[paste0("state", p[1])] <- -1
      est <- coef_of(p[2]) - coef_of(p[1])
      ok <- !is.na(b)
      se <- sqrt(as.numeric(t(cvec[ok]) %*% V[ok, ok, drop = FALSE] %*% cvec[ok]))
      tval <- est / se
      data.frame(region = rname, Deut.Time = tlabel, state1 = p[1], state2 = p[2],
                 estimate = est, se = se,
                 p_value = 2 * stats::pt(-abs(tval), df = stats::df.residual(main)),
                 stringsAsFactors = FALSE)
    }))
    list(summary = summ, pairs = pr)
  }

  res <- list()
  for (rname in names(regions)) {
    r <- regions[[rname]]
    dr <- long[long$Start >= r[1] & long$End <= r[2], , drop = FALSE]
    if (time_model) {
      res[[length(res) + 1]] <- fit_one(dr, rname, "all", TRUE)
    } else {
      for (t in unique(a$Deut.Time)) {
        res[[length(res) + 1]] <- fit_one(dr[dr$Deut.Time == t, , drop = FALSE], rname, t, FALSE)
      }
    }
  }
  summary_df <- do.call(rbind, lapply(res, `[[`, "summary"))
  pairs_df <- do.call(rbind, lapply(res, `[[`, "pairs"))
  if (is.null(summary_df)) stop("no region has two or more peptides present in every state", call. = FALSE)
  summary_df$p_state_adj <- stats::p.adjust(summary_df$p_state, method = p_adjust)
  first <- c("region", "Deut.Time", "n_peptides", "n_obs", "F_state", "df_state", "df_resid",
             "p_state", "p_state_adj")
  summary_df <- summary_df[, c(first, setdiff(colnames(summary_df), first))]
  if (!is.null(pairs_df)) pairs_df$p_adj <- stats::p.adjust(pairs_df$p_value, method = p_adjust)
  rownames(summary_df) <- NULL
  if (!is.null(pairs_df)) rownames(pairs_df) <- NULL
  list(summary = summary_df, pairs = pairs_df)
}

## p-value of the nested-model F test small vs big; NA when not testable
.anova_p <- function(small, big) {
  av <- tryCatch(stats::anova(small, big), error = function(e) NULL)
  if (is.null(av) || nrow(av) < 2 || is.na(av$Df[2]) || av$Df[2] == 0 || av$Res.Df[2] == 0) return(NA_real_)
  av$`Pr(>F)`[2]
}
