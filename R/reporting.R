## Reporting: summary table (Masson et al. 2019), tidy results table,
## provenance (parameters, package version, session).

## attach analysis parameters and versions to a result
.hdx_stamp <- function(x, params) {
  attr(x, "hdxboxer") <- list(
    package_version = as.character(utils::packageVersion("HDXBoxeR")),
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    R_version = R.version.string,
    parameters = params)
  x
}

## write a "<file>_parameters.txt" next to an exported table
.hdx_write_params <- function(file, stamp) {
  pfile <- sub("(\\.[A-Za-z0-9]+)?$", "_parameters.txt", file)
  p <- stamp$parameters
  val <- vapply(p, function(v) {
    s <- trimws(as.character(v))
    if (!is.null(names(v))) s <- paste0(names(v), "=", s)
    paste(s, collapse = ", ")
  }, character(1))
  lines <- c(sprintf("HDXBoxeR %s", stamp$package_version), sprintf("date: %s", stamp$date),
             sprintf("%s", stamp$R_version), "", "parameters:",
             sprintf("  %s: %s", names(val), val), "", "sessionInfo:",
             utils::capture.output(print(utils::sessionInfo())))
  writeLines(lines, pfile)
  invisible(pfile)
}

#' Tidy table of peptide-level results
#'
#' One row per peptide, time point and compared pair of states, with the
#' means, SDs, replicate counts, difference, Welch p-value, adjusted
#' p-value, critical interval and significance call - the package-wide rule
#' (p < `pvalue` AND |difference| > pooled CI). Suitable for supplementary
#' tables and for further analysis outside HDXBoxeR. The analysis
#' parameters, package version and date are attached as attribute
#' `"hdxboxer"`; with `file`, the table is written as csv together with a
#' `<file>_parameters.txt` holding the parameters and `sessionInfo()`.
#'
#' @param filepath input file or table, as for [output_tp()].
#' @param states states to compare. Default: all.
#' @param times time points. Default: all.
#' @param reference reference state; every other state is compared with it
#'   (difference = other - reference). Default: the first state.
#'   `reference = NULL` or `all_pairs = TRUE` compares every pair.
#' @param all_pairs `TRUE` compares every pair of states instead.
#' @param alpha alpha of the pooled critical interval. Default `0.01`.
#' @param pvalue p-value cutoff. Default `0.01`.
#' @param p_adjust multiple-testing correction passed to [stats::p.adjust()]:
#'   `"none"` (default) or e.g. `"BH"`, across the peptides of each pair of
#'   states at each time point.
#' @param percent value scale as in [output_tp()]: `FALSE` (Da, default),
#'   `TRUE`, `"FD"` or `"per_residue"`.
#' @param replicates optional maximum number of replicates per state.
#' @param file optional csv file to write.
#' @return data frame (invisibly when `file` is given).
#' @seealso [hdx_summary_table()], [output_prep()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' res <- hdx_results_table(file_nm)
#' head(res[res$is_sig, ])
#' attr(res, "hdxboxer")$parameters$p_adjust
#' @export
hdx_results_table <- function(filepath, states, times, reference, all_pairs = FALSE,
                              alpha = 0.01, pvalue = 0.01, p_adjust = "none",
                              percent = FALSE, replicates = NULL, file = NULL) {
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  raw <- .hdx_raw(filepath)
  if (missing(states)) states <- unique(arg_df(raw)$Protein.State)
  if (length(states) < 2) stop("need at least two states", call. = FALSE)
  a <- suppressMessages(if (missing(times)) {
    output_tp(raw, states = states, replicates = replicates, percent = percent)
  } else output_tp(raw, states = states, times = times, replicates = replicates, percent = percent))
  if (missing(reference)) reference <- states[1]
  pairs <- if (isTRUE(all_pairs) || is.null(reference)) utils::combn(states, 2, simplify = FALSE) else {
    if (!reference %in% states) stop("reference '", reference, "' not among states", call. = FALSE)
    lapply(setdiff(states, reference), function(s) c(reference, s))
  }

  rows <- list()
  for (t in unique(a$Deut.Time)) {
    sub <- a[a$Deut.Time == t, , drop = FALSE]
    for (pr in pairs) {
      m1 <- .hdx_block_matrix(sub, .hdx_state_cols(pr[1], colnames(sub)))
      m2 <- .hdx_block_matrix(sub, .hdx_state_cols(pr[2], colnames(sub)))
      s1 <- .hdx_row_sd(m1); s2 <- .hdx_row_sd(m2); n1 <- .hdx_row_n(m1); n2 <- .hdx_row_n(m2)
      ci <- CI_2pts(s1, s2, alpha = alpha, n1 = n1, n2 = n2)
      rows[[length(rows) + 1]] <- data.frame(
        Deut.Time = t, Start = sub$Start, End = sub$End, Sequence = sub$Sequence, Charge = sub$Charge,
        state1 = pr[1], state2 = pr[2],
        mean1 = .hdx_row_mean(m1), sd1 = s1, n1 = n1, mean2 = .hdx_row_mean(m2), sd2 = s2, n2 = n2,
        stringsAsFactors = FALSE)
      rows[[length(rows)]]$delta <- rows[[length(rows)]]$mean2 - rows[[length(rows)]]$mean1
      rows[[length(rows)]]$p_value <- .hdx_row_welch(m1, m2)
      rows[[length(rows)]]$ci <- ci
    }
  }
  out <- do.call(rbind, rows)
  out$p_adj <- stats::ave(out$p_value, out$Deut.Time, out$state1, out$state2,
                          FUN = function(p) stats::p.adjust(p, method = p_adjust))
  out$is_sig <- !is.na(out$p_adj) & out$p_adj < pvalue & !is.na(out$delta) & abs(out$delta) > out$ci
  out <- out[, c("Deut.Time", "Start", "End", "Sequence", "Charge", "state1", "state2", "mean1", "sd1", "n1",
                 "mean2", "sd2", "n2", "delta", "p_value", "p_adj", "ci", "is_sig")]
  rownames(out) <- NULL
  out <- .hdx_stamp(out, list(states = states, times = unique(a$Deut.Time),
                              pairs = vapply(pairs, paste, character(1), collapse = " vs "),
                              alpha = alpha, pvalue = pvalue, p_adjust = p_adjust,
                              percent = format(percent), replicates = attr(a, "replicates"),
                              rule = "p_adj < pvalue AND |delta| > pooled CI (Hageman & Weis)"))
  if (!is.null(file)) {
    utils::write.csv(out, file, row.names = FALSE)
    .hdx_write_params(file, attr(out, "hdxboxer"))
    return(invisible(out))
  }
  out
}

#' One-command HTML report of an HDX-MS comparison
#'
#' Renders an HTML report with the analysis parameters, the Masson et al.
#' (2019) summary table ([hdx_summary_table()]), data-quality checks
#' ([peptide_qc()] issues and the coverage map), significant-peptide counts
#' and a fissure-vent plot for every time point and state compared with the
#' reference, and - when given - region ANOVAs ([compare_regions()]) and
#' whole-time-course results ([compare_timecourse()]). The tidy results
#' table ([hdx_results_table()]) is written next to the report as csv,
#' together with its parameters file.
#'
#' Needs the rmarkdown and knitr packages (and pandoc, which comes with
#' RStudio).
#'
#' @param filepath input file or table, as for [output_tp()].
#' @param output_file path of the HTML report. Default `"hdx_report.html"`
#'   in the working directory.
#' @param states,times,reference states, time points and reference state
#'   (default: all states, all common time points, the first state).
#' @param regions optional named list of `c(start, end)` residue ranges for
#'   [compare_regions()].
#' @param alpha,pvalue,p_adjust,percent analysis settings, as in
#'   [hdx_results_table()].
#' @param metadata named list of extra summary-table items, as in
#'   [hdx_summary_table()].
#' @param title report title.
#' @param Dfact D2O fraction of the labeling buffer, for back-exchange.
#'   Default `0.85`.
#' @param protein_range residue range for the coverage, as in
#'   [general_info()].
#' @return the path of the report, invisibly.
#' @examples
#' \dontrun{
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' render_report(file_nm, output_file = file.path(tempdir(), "report.html"),
#'               regions = list(nterm = c(1, 100), cterm = c(101, 225)))
#' }
#' @export
render_report <- function(filepath, output_file = "hdx_report.html", states = NULL, times = NULL,
                          reference = NULL, regions = NULL, alpha = 0.01, pvalue = 0.01,
                          p_adjust = "none", percent = FALSE, metadata = list(),
                          title = "HDX-MS report", Dfact = 0.85, protein_range = NULL) {
  for (pkg in c("rmarkdown", "knitr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("render_report() needs the '", pkg, "' package: install.packages('", pkg, "')", call. = FALSE)
    }
  }
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  template <- system.file("report", "hdx_report.Rmd", package = "HDXBoxeR")
  if (!nzchar(template)) stop("report template not found", call. = FALSE)
  work <- tempfile("hdx_report_")
  dir.create(work)
  file.copy(template, work)
  data_rds <- file.path(work, "data.rds")
  saveRDS(.hdx_raw(filepath), data_rds)
  output_file <- normalizePath(output_file, mustWork = FALSE)
  results_file <- sub("\\.html?$", "", output_file)
  results_file <- paste0(results_file, "_results.csv")
  rmarkdown::render(file.path(work, "hdx_report.Rmd"), output_file = basename(output_file),
                    output_dir = dirname(output_file), quiet = TRUE, envir = new.env(),
                    intermediates_dir = work,
                    params = list(data = data_rds, title = title, states = states, times = times,
                                  reference = reference, regions = regions, alpha = alpha,
                                  pvalue = pvalue, p_adjust = p_adjust, percent = percent,
                                  metadata = metadata, results_file = results_file,
                                  Dfact = Dfact, protein_range = protein_range))
  message("report written to ", output_file)
  invisible(output_file)
}

#' Summary table for HDX-MS reporting (Masson et al. 2019)
#'
#' Builds the summary table recommended by Masson et al. (2019, Nat
#' Methods 16:595) for publishing HDX-MS data: one row per item, one column
#' per protein state. Items computed from the data: time course, controls,
#' replicates, number of peptides, sequence coverage, average peptide
#' length, redundancy, back-exchange (mean and range), repeatability
#' (pooled replicate SD, Da) and the significance threshold: the critical
#' interval used when each state is compared with the first state (the
#' same values as in [hdx_results_table()]). Items that only you know
#' (labeling buffer, pH, temperature, quench, protease, LC/MS, software
#' ...) are added from `metadata`.
#'
#' @param filepath input file or table, as for [output_tp()].
#' @param metadata named list of extra items, e.g.
#'   `list("HDX reaction details" = "20 mM HEPES pD 7.4, 90% D2O, 20 C")`.
#'   A single value applies to every state; a vector gives one value per
#'   state. They are listed first.
#' @param replicates optional maximum number of replicates per state.
#'   Default `NULL` uses every replicate.
#' @param Dfact D2O fraction of the labeling buffer for back-exchange.
#'   Default `0.85`.
#' @param alpha alpha of the critical interval. Default `0.01`.
#' @param pvalue p-value cutoff, reported with the threshold. Default `0.01`.
#' @param protein_range residue range for the coverage, as in
#'   [general_info()]. Default: first to last covered residue.
#' @param file optional output file: `.csv` or `.md` (markdown table).
#' @param states optional states to include (the first one is the
#'   reference). Default: every state.
#' @return data frame with columns `Item` and one per state (invisibly
#'   when `file` is given).
#' @seealso [all_summary()], [output_prep()], [hdx_results_table()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' hdx_summary_table(file_nm, metadata = list(
#'   "HDX reaction details" = "PBS pD 7.4, 85% D2O, 22 C",
#'   "Data analysis software" = "HDExaminer 3; HDXBoxeR"))
#' @export
hdx_summary_table <- function(filepath, metadata = list(), replicates = NULL, Dfact = 0.85,
                              alpha = 0.01, pvalue = 0.01, protein_range = NULL, file = NULL,
                              states = NULL) {
  raw <- .hdx_raw(filepath)
  if (!is.null(states)) {
    miss <- setdiff(states, raw$Protein.State)
    if (length(miss)) stop("state(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
    raw <- raw[raw$Protein.State %in% states, , drop = FALSE]
    ## the first state given is the reference: put it first in the data
    raw <- raw[order(match(raw$Protein.State, states)), , drop = FALSE]
  }
  s <- suppressMessages(all_summary(raw, replicates = replicates, Dfact = Dfact, alpha = alpha,
                                    protein_range = protein_range))
  states <- s$Protein.State
  ## significance threshold: the critical interval actually used for each
  ## state vs the first (reference) state, across time points (Da)
  res <- suppressMessages(hdx_results_table(raw, states = states, alpha = alpha, pvalue = pvalue,
                                            replicates = replicates))
  thr <- vapply(states, function(st) {
    if (st == states[1]) return("reference state")
    ci <- unique(round(res$ci[res$state2 == st], 2))
    rng <- if (length(ci) == 1) sprintf("%.2f", ci) else sprintf("%.2f-%.2f", min(ci), max(ci))
    sprintf("|difference| > %s Da (critical interval, alpha = %s) and Welch p < %s", rng, alpha, pvalue)
  }, character(1))
  bx <- ifelse(is.na(s$`<BackHX>`), "not available",
               sprintf("%s (%s - %s)", s$`<BackHX>`, s$`BackHX Range From`, s$To))
  pr <- attr(general_info(raw, protein_range = protein_range), "protein_range")
  ctrl <- vapply(states, function(st) {
    tt <- unique(raw$Deut.Time[raw$Protein.State == st])
    paste(c(if ("0s" %in% tt) "undeuterated (0s)", if ("FD" %in% tt) "fully deuterated (FD)"),
          collapse = "; ")
  }, character(1))
  computed <- list(
    "HDX time course" = s$Timepoints,
    "HDX controls" = ifelse(ctrl == "", "none", ctrl),
    "Number of replicates" = s$`# replicates`,
    "Number of peptides" = s$`# peptides`,
    "Sequence coverage (%)" = sprintf("%s (residues %s-%s)", s$`peptide coverage %`, pr[1], pr[2]),
    "Average peptide length" = s$`<Peptide Length>`,
    "Redundancy" = s$`<Redundancy>`,
    "Back-exchange (%): mean (range)" = bx,
    "Repeatability (pooled SD, Da)" = s$st.dev,
    "Significant differences in HDX" = thr
  )
  items <- c(metadata, computed)
  rows <- lapply(names(items), function(nm) {
    v <- as.character(items[[nm]])
    if (length(v) == 1) v <- rep(v, length(states))
    if (length(v) != length(states)) stop("metadata item '", nm, "' must have 1 or ", length(states), " values", call. = FALSE)
    stats::setNames(as.list(c(nm, v)), c("Item", states))
  })
  out <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
  colnames(out) <- c("Item", states)
  rownames(out) <- NULL
  out <- .hdx_stamp(out, list(replicates = if (is.null(replicates)) "all" else replicates,
                              Dfact = Dfact, alpha = alpha, pvalue = pvalue,
                              protein_range = pr))
  if (!is.null(file)) {
    if (grepl("\\.md$", file, ignore.case = TRUE)) {
      md <- c(paste("|", paste(colnames(out), collapse = " | "), "|"),
              paste("|", paste(rep("---", ncol(out)), collapse = " | "), "|"),
              apply(out, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
      writeLines(md, file)
    } else utils::write.csv(out, file, row.names = FALSE)
    return(invisible(out))
  }
  out
}
