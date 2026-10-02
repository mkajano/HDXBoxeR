#' Compare peptide-level HDX uptake across three or more protein states
#'
#' Every other significance function in this package (`plot_fissure_vent()`,
#' `diff_peptides_prep()`, `robot_indexes()`) compares exactly two states at
#' a time - three or more states means looping pairwise yourself, with no
#' correction for now running many more simultaneous tests. This function
#' does two things instead: (1) an omnibus Kruskal-Wallis test per peptide
#' across *all* requested states (does this peptide differ across any of
#' them at all) - reported, not used to gate significance (see Details) -
#' and (2) the package's usual Hageman & Weis hybrid test (pooled
#' `CI_2pts()` + per-pair Welch p-value) for every pair of states (or every
#' `reference` vs. the rest, if `reference` is given), with p-values
#' adjusted across every peptide-by-pair test at this timepoint together.
#'
#' @section Which replicate SD feeds the pooled CI:
#' `plot_fissure_vent()`'s default `sd_source = "CI"` derives each state's
#' SD from a two-sample CI against `df`'s *first* state - a designated
#' reference baked into that formula, harmless when every comparison is
#' against one fixed reference, but not well-defined for an arbitrary pair
#' among 3+ states with no single reference. This function instead always
#' uses each state's own raw replicate SD directly (`plot_fissure_vent()`'s
#' `sd_source != "CI"` alternative), which treats every state symmetrically
#' and needs no reference concept - the correct generalization here, not
#' merely an implementation convenience. One consequence: `ci` values from
#' this function will not numerically match `plot_fissure_vent()`'s default
#' for the same pair (they will match `plot_fissure_vent(sd_source = "raw SD")`,
#' i.e. any `sd_source` value other than `"CI"`).
#'
#' @section Available-case, not complete-case, matching:
#' `output_tp()`'s wide-format merge is an inner join across every state
#' passed to it - a peptide missing from *any one* of N states is dropped
#' entirely, so calling `output_tp()` once with all N states would silently
#' shrink the comparable peptide set combinatorially as N grows (a peptide
#' identified in 90% of runs per state survives in only 0.9^N of cases).
#' This function instead calls `output_tp()` **separately per pair**, so a
#' peptide only needs to be present in the *two* states of that specific
#' pair to be tested there - "available-case" matching. The one place
#' complete-case is unavoidable is the omnibus test itself: Kruskal-Wallis
#' cannot test a peptide against a group it has no data for, so the omnibus
#' p-value is computed from the peptide set present in *all* requested
#' states (via one additional `output_tp()` call with every state), and
#' peptides outside that set simply get `NA` in `p_omnibus`/`p_omnibus_adj`.
#'
#' @section Why the omnibus test is advisory, not a gate:
#' With the small replicate counts typical here (often `replicates = 3`),
#' Kruskal-Wallis's exact null distribution is extremely coarse (with 3
#' states of 3 replicates each, the smallest attainable p-value is
#' 6/1680 ~ 0.0036) and its chi-squared approximation is unreliable at this
#' size. Using it as a hard gate stacks a second, statistically shaky
#' filter on top of the BH-adjusted pairwise p-value with no clear combined
#' error-rate guarantee. `is_sig` therefore depends only on the per-pair
#' adjusted p-value and the pooled CI - the direct N-state generalization of
#' the shipped two-state criterion - and `p_omnibus`/`p_omnibus_adj` are
#' reported alongside for context, not enforced.
#'
#' @param filepath filepath to the input file (same format as
#'   [output_tp()]/[nm_states()] expect) - this function calls
#'   [output_tp()] internally, once per pair plus once for all states
#'   together, rather than taking an already-built wide table, precisely
#'   because that wide format has already discarded the availability
#'   information available-case matching needs (see Details).
#' @param times a single `Deut.Time` value to compare at.
#' @param states character vector of three or more states to compare.
#'   Default `NULL` uses every state found in `filepath` ([nm_states()]).
#' @param by matching rule for identifying "the same peptide" across
#'   states. Currently only `"sequence"` (Sequence+Start+End identity,
#'   [output_tp()]'s `seq_match = TRUE`) is implemented. A `"position"`
#'   mode (Start+End only, ignoring Sequence - needed so a peptide covering
#'   a point mutation isn't dropped just because its sequence differs
#'   between states) is planned but not yet implemented; passing it errors
#'   rather than silently falling back.
#' @param reference optional state name. If given, only `reference` vs.
#'   each other state is tested (N-1 pairs, Dunnett-style), instead of
#'   every pair (N(N-1)/2) - both fewer simultaneous tests and usually the
#'   actually-relevant contrast for a wild-type/apo vs. several perturbations
#'   design. Default `NULL` tests every pair.
#' @param regions optional named list of `c(start, end)` residue ranges
#'   (e.g. `list(NTR = c(1, 85), ACD = c(86, 165))`). When supplied, only
#'   peptides fully contained in a region (`Start >= start & End <= end`)
#'   are compared, run separately per region, and the result gets a
#'   `region` column. Default `NULL` compares every peptide with no
#'   filtering (no `region` column in the output).
#' @param replicates number of replicates per state. Default `3`.
#' @param alpha alpha for the pooled critical interval (`CI_2pts()`),
#'   applied per pair. Default `0.01`.
#' @param pvalue significance cutoff applied to the (BH-adjusted) per-pair
#'   p-value. Default `0.01`.
#' @param p_adjust_method passed to [stats::p.adjust()]; applied separately
#'   to the vector of per-peptide omnibus p-values and the vector of
#'   per-peptide-by-pair p-values (the latter across every pair together,
#'   not within each pair on its own). Default `"BH"`.
#' @param percent passed to [output_tp()] - `FALSE` (default) compares
#'   uptake in Daltons, `TRUE` compares percent deuteration.
#' @return data.frame with one row per peptide-by-state-pair (by region, if
#'   `regions` is supplied): `peptide_id`, `Start`, `End`, `Sequence`,
#'   `state1`, `state2`, `delta` (`state2 - state1`), `p_value` (raw
#'   per-pair Welch p-value), `ci` (pooled CI for that pair), `p_adj`
#'   (BH-adjusted per-pair p-value), `p_omnibus`/`p_omnibus_adj` (`NA` for
#'   peptides absent from at least one requested state), and `is_sig`.
#' @seealso [plot_fissure_vent()], [diff_peptides_prep()], [CI_2pts()],
#'   [compare_states_unmatched()]
#' @examples
#' \donttest{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' # the shipped example only has two states; compare_states_tp() needs 3+.
#' # See vignette("anova") for a worked example that
#' # constructs a synthetic third state to demonstrate this function.
#' }
#' @export
compare_states_tp <- function(filepath, times, states = NULL,
                               by = c("sequence", "position"),
                               reference = NULL, regions = NULL,
                               replicates = 3, alpha = 0.01, pvalue = 0.01,
                               p_adjust_method = "BH", percent = FALSE) {
  by <- match.arg(by)
  if (by == "position") {
    stop("by = \"position\" is not yet implemented - only by = \"sequence\" (the default) works today.")
  }

  if (is.null(states)) states <- nm_states(filepath)
  if (length(states) < 3) {
    stop("compare_states_tp() is for three or more states - for exactly two, use plot_fissure_vent()/diff_peptides_prep().")
  }
  if (!is.null(reference) && !(reference %in% states)) {
    stop(sprintf("reference '%s' not among states: %s", reference, paste(states, collapse = ", ")))
  }

  pairs <- if (!is.null(reference)) {
    lapply(setdiff(states, reference), function(s) c(reference, s))
  } else {
    utils::combn(states, 2, simplify = FALSE)
  }

  in_region <- function(d, region) {
    if (is.null(region)) return(rep(TRUE, nrow(d)))
    d$Start >= region[1] & d$End <= region[2]
  }
  region_list <- if (is.null(regions)) list(all = NULL) else regions

  ## output_tp() is called once per pair below (plus once for the omnibus) -
  ## deliberately (see "Available-case, not complete-case, matching" above),
  ## not something to collapse into a single call. What CAN be shared is the
  ## raw file parse itself: arg_df(filepath) is a pure function of the file
  ## and doesn't depend on which states/pair are being extracted, so reading
  ## and parsing it once here and passing it into every output_tp() call
  ## (its .arg_df argument) turns N(N-1)/2 + 1 full re-reads of the same
  ## file into exactly one, with identical results - a real cost for
  ## several states or a large file, and pure overhead otherwise.
  raw_arg_df <- arg_df(filepath)

  ## omnibus: needs the peptide present in every requested state at once -
  ## the one place complete-case is statistically unavoidable (see Details)
  omni_df <- output_tp(filepath, replicates = replicates, states = states, times = times,
                        seq_match = TRUE, percent = percent, .arg_df = raw_arg_df)
  omni_df <- omni_df[omni_df$Deut.Time == times, ]
  peptide_id_omni <- paste(omni_df$Start, omni_df$End, omni_df$Sequence, sep = "_")

  p_omnibus_raw <- numeric(0)
  if (nrow(omni_df) > 0) {
    value_cols_omni <- colnames(omni_df)[7:ncol(omni_df)]
    state_of_col_omni <- sub("_(X\\.\\.Deut|Deut\\.\\.)_[0-9]+$", "", value_cols_omni)
    cols_for_omni <- function(s) value_cols_omni[state_of_col_omni == .state_column_prefix(s)]
    p_omnibus_raw <- vapply(seq_len(nrow(omni_df)), function(r) {
      vals <- unlist(lapply(states, function(s) as.numeric(omni_df[r, cols_for_omni(s)])))
      grp  <- rep(states, each = replicates)
      ok <- !is.na(vals)
      if (length(unique(grp[ok])) < length(states) || length(unique(vals[ok])) < 2) return(NA_real_)
      tryCatch(stats::kruskal.test(vals[ok], factor(grp[ok]))$p.value, error = function(e) NA_real_)
    }, numeric(1))
  }

  ## per-pair: same hybrid test as plot_fissure_vent()/diff_peptides_prep(),
  ## generalized to every requested pair. Computed fully self-contained
  ## (each state's own replicate SD directly, pooled via CI_2pts()) rather
  ## than via sd_ci_timepoint(), which assumes the data.frame's own first
  ## state-block is the reference for every column it returns - true for
  ## the shipped two-state functions, not safe to assume for an arbitrary
  ## pair here.
  build_pair <- function(pr) {
    s1 <- pr[1]; s2 <- pr[2]
    pair_df <- output_tp(filepath, replicates = replicates, states = c(s1, s2), times = times,
                          seq_match = TRUE, percent = percent, .arg_df = raw_arg_df)
    pair_df <- pair_df[pair_df$Deut.Time == times, ]
    if (nrow(pair_df) == 0) return(NULL)

    value_cols <- colnames(pair_df)[7:ncol(pair_df)]
    state_of_col <- sub("_(X\\.\\.Deut|Deut\\.\\.)_[0-9]+$", "", value_cols)
    raw1 <- value_cols[state_of_col == .state_column_prefix(s1)]
    raw2 <- value_cols[state_of_col == .state_column_prefix(s2)]
    if (length(raw1) == 0 || length(raw2) == 0) {
      stop(sprintf("could not match state '%s' or '%s' to any column among: %s - this should not happen, please file a bug report",
                    s1, s2, paste(value_cols, collapse = ", ")))
    }

    sd1 <- apply(as.matrix(pair_df[, raw1]), 1, stats::sd, na.rm = TRUE)
    sd2 <- apply(as.matrix(pair_df[, raw2]), 1, stats::sd, na.rm = TRUE)
    ci  <- CI_2pts(stats::na.omit(sd1), stats::na.omit(sd2), replicates, alpha)

    test <- lapply(seq_len(nrow(pair_df)), function(r) {
      x1 <- as.numeric(pair_df[r, raw1]); x2 <- as.numeric(pair_df[r, raw2])
      tt <- tryCatch(stats::t.test(x1, x2), error = function(e) NULL)
      if (is.null(tt) || is.na(tt$p.value)) return(c(delta = NA_real_, p_value = NA_real_))
      c(delta = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE), p_value = tt$p.value)
    })
    test <- do.call(rbind, test)

    peptide_id <- paste(pair_df$Start, pair_df$End, pair_df$Sequence, sep = "_")
    omni_match <- match(peptide_id, peptide_id_omni)

    data.frame(
      peptide_id = peptide_id, Start = pair_df$Start, End = pair_df$End, Sequence = pair_df$Sequence,
      state1 = s1, state2 = s2,
      delta = test[, "delta"], p_value = test[, "p_value"], ci = ci,
      p_omnibus = if (length(p_omnibus_raw)) p_omnibus_raw[omni_match] else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  ## build_pair() (a full output_tp() call plus a t-test per peptide) does
  ## not depend on `region` at all - only the filter below does - so it is
  ## computed exactly once per pair here and reused across every region,
  ## rather than once per region-by-pair combination.
  pair_results <- lapply(pairs, build_pair)

  out_all <- lapply(names(region_list), function(rname) {
    region <- region_list[[rname]]
    pair_out <- lapply(pair_results, function(d) {
      if (is.null(d)) return(NULL)
      d[in_region(d, region), ]
    })
    d <- do.call(rbind, pair_out)
    if (!is.null(d)) d$region <- rname
    d
  })
  out <- do.call(rbind, out_all)
  if (is.null(out) || nrow(out) == 0) {
    stop(sprintf("no peptides found for any pair at Deut.Time == '%s'", times))
  }

  ## BH over every peptide x pair test together, THEN gate on the pooled
  ## CI - never subset to CI-passing peptides before adjusting, or the
  ## realized FDR is no longer what p.adjust() promises.
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust_method)
  uniq <- !duplicated(out$peptide_id)
  omni_adj_uniq <- stats::p.adjust(out$p_omnibus[uniq], method = p_adjust_method)
  out$p_omnibus_adj <- omni_adj_uniq[match(out$peptide_id, out$peptide_id[uniq])]

  out$is_sig <- !is.na(out$p_adj) & (abs(out$delta) > out$ci) & (out$p_adj < pvalue)
  if (is.null(regions)) out$region <- NULL
  rownames(out) <- NULL
  out
}
