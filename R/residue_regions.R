#' Cluster significant peptides into contiguous destabilized/protected regions
#'
#' Takes a peptide-level comparison between two states at one timepoint and
#' derives contiguous residue regions that are consistently destabilized
#' (more exposed) or protected (less exposed), rather than reporting a flat
#' list of significant peptides. A peptide-level significance call is
#' length-biased - a long peptide can clear a flat critical-interval cutoff
#' from a diffuse, unremarkable per-residue effect just by accumulating a
#' large summed Da difference - so this function gates on FDR-corrected
#' per-peptide p-values (scale-invariant to peptide length) instead, ranks
#' by per-residue effect size (delta divided by peptide length), and only
#' lets a residue into a region once at least `min_peptides` independent
#' peptides of the *same direction* cover it and clear the FDR cutoff.
#' Requiring multiple corroborating peptides guards against a single long
#' peptide alone driving a region's inclusion.
#'
#' @param df output of [output_tp()] (or `output_tp_proc()`) called with
#'   `percent = FALSE`, i.e. deuteration uptake in Daltons - this function
#'   works in Da so that `delta` is additive across residues, same
#'   requirement as [residue_deconvolution()].
#' @param state1 First protein state, matched the same way as in
#'   [plot_fissure_vent()] (substring match against the raw replicate
#'   column names - does not need to be the exact column name).
#' @param state2 Second protein state, same matching rule.
#' @param times a single Deut.Time value to derive regions for (e.g.
#'   `"1800.00s"`). Call this function once per timepoint of interest and
#'   `rbind()` the results if you want regions across multiple timepoints.
#' @param replicates number of replicates in sample. Default set to 3.
#' @param fdr_cutoff BH-adjusted p-value cutoff for a peptide to count as
#'   significant. Default `0.05`. Unlike [plot_fissure_vent()]'s per-peptide
#'   p-value (a raw, uncorrected Welch p-value), this is FDR-corrected
#'   across every peptide tested in this comparison, since region-finding
#'   involves many more simultaneous per-peptide tests than eyeballing one
#'   plot.
#' @param min_peptides minimum number of independent, same-direction,
#'   FDR-significant peptides that must cover a residue before that residue
#'   is included in a region. Default `2`.
#' @param gap_tolerance maximum number of consecutive non-covered residues
#'   allowed inside one region before it is split into two. Default `2`.
#' @param alpha alpha used only for the diagnostic pooled Hageman & Weis
#'   critical interval reported alongside each region (via
#'   [sd_ci_timepoint()]/[CI_2pts()], the same pooled threshold
#'   [plot_fissure_vent()] draws as its grey band) - not used as the
#'   region-inclusion gate itself. Default `0.01`.
#' @return data.frame with one row per region: `direction`
#'   ("destabilized"/"protected"), `region_start`, `region_end`,
#'   `n_residues`, `n_peptides` (distinct FDR-significant, same-direction
#'   peptides overlapping the region), `mean_per_residue_delta` (Da per
#'   residue), and `pooled_CI_Da` (the diagnostic Hageman & Weis pooled
#'   threshold for this comparison, repeated on every row for reference).
#' @seealso [plot_fissure_vent()], [residue_deconvolution()],
#'   [residue_level_average()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm, percent = FALSE)
#' residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' @export
residue_regions <- function(df, state1, state2, times, replicates = 3,
                             fdr_cutoff = 0.05, min_peptides = 2,
                             gap_tolerance = 2, alpha = 0.01) {

  sub_df <- df[df$Deut.Time == times, ]
  if (nrow(sub_df) == 0) {
    stop(sprintf("no rows with Deut.Time == '%s'", times))
  }

  resolve_raw_cols <- function(state, cols, replicates) {
    state <- .state_column_prefix(state)
    hits <- cols[startsWith(cols, state)]
    if (length(hits) != replicates) {
      stop(sprintf("Expected %d raw replicate columns matching state '%s', found %d: %s",
                    replicates, state, length(hits), paste(hits, collapse = ", ")))
    }
    hits
  }
  raw_cols1 <- resolve_raw_cols(state1, colnames(sub_df), replicates)
  raw_cols2 <- resolve_raw_cols(state2, colnames(sub_df), replicates)

  # per-peptide delta (state2 - state1, Da) and Welch p-value, same
  # per-peptide test as plot_fissure_vent()'s significance half - computed
  # independently here (not shared code) since this function only needs
  # the p-value and point delta, not a calibrated confidence interval.
  test <- lapply(seq_len(nrow(sub_df)), function(r) {
    x1 <- as.numeric(sub_df[r, raw_cols1]); x2 <- as.numeric(sub_df[r, raw_cols2])
    tt <- tryCatch(t.test(x1, x2), error = function(e) NULL)
    if (is.null(tt) || is.na(tt$p.value)) {
      return(c(delta = NA_real_, p_value = NA_real_))
    }
    c(delta = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE), p_value = tt$p.value)
  })
  test <- do.call(rbind, test)

  pep <- data.frame(sub_df[, c("Start", "End")], delta = test[, "delta"], p_value = test[, "p_value"])
  pep$length <- pep$End - pep$Start + 1
  pep$per_residue_delta <- pep$delta / pep$length
  pep$p_adj <- stats::p.adjust(pep$p_value, method = "BH")
  pep$fdr_sig <- !is.na(pep$p_adj) & pep$p_adj < fdr_cutoff

  max_pos <- max(pep$End)

  # for each residue, how many independent FDR-significant peptides of the
  # given direction cover it, and their mean per-residue delta
  residue_hit_counts <- function(direction) {
    counts <- integer(max_pos)
    sum_prd <- numeric(max_pos)
    hit <- pep[pep$fdr_sig & (if (direction == "destabilized") pep$per_residue_delta > 0
                               else pep$per_residue_delta < 0), ]
    for (i in seq_len(nrow(hit))) {
      rng <- hit$Start[i]:hit$End[i]
      rng <- rng[rng >= 1 & rng <= max_pos]
      counts[rng] <- counts[rng] + 1
      sum_prd[rng] <- sum_prd[rng] + hit$per_residue_delta[i]
    }
    list(counts = counts, mean_prd = ifelse(counts > 0, sum_prd / counts, NA_real_))
  }

  cluster_regions <- function(confirmed_positions) {
    if (length(confirmed_positions) == 0) return(list())
    pos <- sort(confirmed_positions)
    breaks <- which(diff(pos) > gap_tolerance + 1)
    starts <- c(1, breaks + 1)
    ends <- c(breaks, length(pos))
    lapply(seq_along(starts), function(i) c(pos[starts[i]], pos[ends[i]]))
  }

  sd1 <- sd_ci_timepoint(sub_df, replicates, alpha)
  sd_prefixes <- c("sd_ci_", "sd_")
  resolve_col <- function(state, cols, known_prefixes) {
    state <- .state_column_prefix(state)
    root <- cols
    for (p in known_prefixes) root <- sub(paste0("^", p), "", root)
    cols[startsWith(root, state)][1]
  }
  sd_state1 <- resolve_col(state1, colnames(sd1), sd_prefixes)
  sd_state2 <- resolve_col(state2, colnames(sd1), sd_prefixes)
  pooled_CI <- CI_2pts(stats::na.omit(sd1[[sd_state1]]), stats::na.omit(sd1[[sd_state2]]), replicates, alpha)

  out <- list()
  for (direction in c("destabilized", "protected")) {
    hc <- residue_hit_counts(direction)
    confirmed <- which(hc$counts >= min_peptides)
    regions <- cluster_regions(confirmed)

    for (r in regions) {
      lo <- r[1]; hi <- r[2]
      rng <- lo:hi
      sign_ok <- if (direction == "destabilized") `>` else `<`
      covering <- pep[pep$End >= lo & pep$Start <= hi & pep$fdr_sig & sign_ok(pep$per_residue_delta, 0), ]
      out[[length(out) + 1]] <- data.frame(
        direction = direction, region_start = lo, region_end = hi,
        n_residues = length(rng), n_peptides = nrow(covering),
        mean_per_residue_delta = mean(hc$mean_prd[rng], na.rm = TRUE),
        pooled_CI_Da = pooled_CI
      )
    }
  }

  if (length(out) == 0) {
    return(data.frame(direction = character(0), region_start = integer(0), region_end = integer(0),
                       n_residues = integer(0), n_peptides = integer(0),
                       mean_per_residue_delta = numeric(0), pooled_CI_Da = numeric(0)))
  }
  result <- do.call(rbind, out)
  result[order(result$direction, result$region_start), ]
}

#' Roll up derived regions against named residue subsets
#'
#' A companion to [residue_regions()]: that function always scans every
#' residue covered by any tested peptide - the unrestricted, data-driven
#' "total" picture, deliberately not pre-restricted to any subset you might
#' already suspect matters. This function takes that unrestricted result
#' and a named list of residue subsets (e.g. structural domains) and
#' reports, per subset, how much of it the data-driven regions actually
#' cover - a cross-check for a structural hypothesis against the
#' unconstrained result, not a replacement for it. Call [residue_regions()]
#' once (it does the expensive per-peptide testing and clustering), then
#' this function as many times as you like against different subset
#' definitions - it's a cheap rollup over an already-computed result.
#'
#' @param region_df the data.frame returned by [residue_regions()].
#' @param regions named list of `c(start, end)` residue ranges to check
#'   against `region_df`, e.g. `list(NTR = c(1, 85), ACD = c(86, 165))`.
#' @return data.frame with one row per named subset: `subset`,
#'   `n_residues` (subset size), `n_residues_significant` (residues in the
#'   subset covered by any region in `region_df`), `frac_significant`,
#'   `n_regions_overlapping`, and `directions` (comma-separated set of
#'   `direction` values among the overlapping regions, `""` if none).
#' @seealso [residue_regions()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' a <- output_tp(file_nm, percent = FALSE)
#' rr <- residue_regions(a, state1 = "Unbound", state2 = "bound", times = "1800.00s")
#' residue_regions_subset_summary(rr, regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' @export
residue_regions_subset_summary <- function(region_df, regions) {
  if (!is.data.frame(region_df) || !all(c("region_start", "region_end", "direction") %in% colnames(region_df))) {
    stop("region_df must be the data.frame returned by residue_regions().")
  }
  if (missing(regions) || is.null(regions) || length(regions) == 0) {
    stop("regions is required - a named list of c(start, end) residue ranges to check.")
  }
  out <- do.call(rbind, lapply(names(regions), function(rname) {
    rng <- regions[[rname]]
    res_in_subset <- rng[1]:rng[2]
    overlap <- region_df[region_df$region_end >= rng[1] & region_df$region_start <= rng[2], ]
    sig_positions <- if (nrow(overlap) > 0) {
      unique(unlist(lapply(seq_len(nrow(overlap)), function(i) overlap$region_start[i]:overlap$region_end[i])))
    } else integer(0)
    n_sig <- sum(res_in_subset %in% sig_positions)
    data.frame(subset = rname, n_residues = length(res_in_subset), n_residues_significant = n_sig,
               frac_significant = n_sig / length(res_in_subset), n_regions_overlapping = nrow(overlap),
               directions = paste(sort(unique(overlap$direction)), collapse = ","),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
