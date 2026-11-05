## Internal helpers for the wide tables made by output_tp()/output_tc().
##
## Columns 1-6 are peptide descriptors; columns 7+ are replicate values named
## "<root>_<k>", where <root> is "<state>_X..Deut" / "<state>_Deut.." for
## output_tp() and "t<time>_X..Deut" / "t<time>_Deut.." for output_tc(), and
## <k> is the replicate number. Blocks (one per state or time) are found from
## these names, so states may have different numbers of replicates and
## missing replicates may be NA. Tables whose value columns do not follow the
## naming scheme fall back to fixed-width blocks of `replicates` columns.

.hdx_blocks <- function(df, replicates = NULL) {
  value_idx <- seq.int(7, ncol(df))
  nms <- colnames(df)[value_idx]
  if (length(nms) > 0 && all(grepl("_[0-9]+$", nms))) {
    roots <- sub("_[0-9]+$", "", nms)
    roots_u <- unique(roots)
    if (!is.unsorted(match(roots, roots_u))) {
      return(stats::setNames(lapply(roots_u, function(r) value_idx[roots == r]), roots_u))
    }
  }
  if (is.null(replicates)) replicates <- 3
  replicates <- replicates[1]
  nb_sets <- length(value_idx) %/% replicates
  if (nb_sets * replicates != length(value_idx)) {
    stop(sprintf("cannot split %d value columns into blocks of %d replicates",
                 length(value_idx), replicates), call. = FALSE)
  }
  idx <- lapply(seq_len(nb_sets), function(j) value_idx[(j - 1) * replicates + seq_len(replicates)])
  roots <- sub("_[0-9]+$", "", nms[(seq_len(nb_sets) - 1) * replicates + 1])
  stats::setNames(idx, roots)
}

## numeric matrix of one block
.hdx_block_matrix <- function(df, cols) {
  m <- as.matrix(df[, cols, drop = FALSE])
  storage.mode(m) <- "double"
  m
}

.hdx_row_sd <- function(m) {
  n <- rowSums(!is.na(m))
  mu <- rowMeans(m, na.rm = TRUE)
  ss <- rowSums((m - mu)^2, na.rm = TRUE)
  out <- sqrt(ss / (n - 1))
  out[n < 2] <- NA_real_
  out
}

.hdx_row_mean <- function(m) {
  out <- rowMeans(m, na.rm = TRUE)
  out[is.nan(out)] <- NA_real_
  out
}

## Welch t-test per row; NA when either side has < 2 values or the test is
## undefined (e.g. all values identical)
.hdx_row_welch <- function(m1, m2, conf.level = 0.95) {
  vapply(seq_len(nrow(m1)), function(i) {
    x1 <- m1[i, ][!is.na(m1[i, ])]; x2 <- m2[i, ][!is.na(m2[i, ])]
    if (length(x1) < 2 || length(x2) < 2) return(NA_real_)
    p <- tryCatch(stats::t.test(x1, x2, conf.level = conf.level)$p.value, error = function(e) NA_real_)
    if (is.nan(p)) NA_real_ else p
  }, numeric(1))
}

## the _proc functions use a percent table and an uptake table row by row:
## they must describe the same peptides in the same order
.hdx_check_aligned <- function(input_proc, input_up) {
  key <- function(d) paste(d$Deut.Time, d$Start, d$End, d$Sequence, d$Charge)
  if (nrow(input_proc) != nrow(input_up) || !identical(key(input_proc), key(input_up))) {
    stop("input_proc and input_up must contain the same peptides in the same order: ",
         "make both with output_tp() from the same data, states and times ",
         "(percent = TRUE and percent = FALSE).", call. = FALSE)
  }
  invisible(TRUE)
}

## warning for a renamed argument (old name still accepted)
.hdx_deprecated <- function(old, new) {
  warning(sprintf("argument '%s' is deprecated; use '%s' instead.", old, new), call. = FALSE)
}

## Multiple-testing adjustment of a p-value table (columns 7+ hold p-values).
## Family, everywhere in the package: the peptides of one pair of states at
## one time point. For pv_timepoint() tables each column is one pair and
## rows are split by Deut.Time; for pv_timecourse() tables each column is
## one time point of one pair. Column 7 of a pv_timepoint() table is the
## reference state compared with itself and is left out.
.hdx_adjust_pv <- function(pv, p_adjust = "none", skip_first = TRUE) {
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  if (p_adjust == "none") return(pv)
  cols <- seq.int(if (skip_first) 8 else 7, ncol(pv))
  if (length(cols) == 0 || cols[1] > ncol(pv)) return(pv)
  groups <- if ("Deut.Time" %in% colnames(pv)) split(seq_len(nrow(pv)), pv$Deut.Time) else list(seq_len(nrow(pv)))
  for (g in groups) for (j in cols) {
    pv[g, j] <- stats::p.adjust(pv[g, j], method = p_adjust)
  }
  pv
}

## per-row count of non-missing replicates in a block
.hdx_row_n <- function(m) rowSums(!is.na(m))

## replicate count per block, as stored by output_tp()/output_tc() or derived
.hdx_replicates <- function(df, replicates = NULL) {
  r <- attr(df, "replicates", exact = TRUE)
  if (!is.null(r)) return(r)
  b <- .hdx_blocks(df, replicates)
  vapply(b, length, integer(1))
}

## SD tables from sd_timepoint() carry attribute "n" (per-peptide replicate
## counts, one column per block, rownames matching the table). Class
## "hdx_sd" keeps that attribute aligned when the table is subset.

#' @export
`[.hdx_sd` <- function(x, ...) {
  n <- attr(x, "n", exact = TRUE)
  out <- NextMethod()
  if (!is.data.frame(out) || is.null(n)) return(out)
  n <- n[match(rownames(out), rownames(x)), , drop = FALSE]
  rownames(n) <- rownames(out)
  roots <- sub("^sd_", "", colnames(out))
  attr(out, "n") <- n[, intersect(colnames(n), roots), drop = FALSE]
  class(out) <- class(x)
  out
}

## per-peptide replicate counts behind SD column `col` of an SD table (rows
## optional); NULL when the table carries no counts
.hdx_n_for <- function(sd_df, col, rows = NULL) {
  n <- attr(sd_df, "n", exact = TRUE)
  if (is.null(n)) return(NULL)
  root <- sub("^sd_", "", col)
  if (!root %in% colnames(n)) return(NULL)
  v <- n[[root]]
  if (!is.null(rows)) v <- v[rows]
  v
}

## pooled SD of one state's per-peptide SDs, weighted by degrees of freedom
## (n_i - 1); with equal n this is sqrt(mean(s^2)), the Hageman & Weis form
.hdx_pooled_sd <- function(s, n) {
  ok <- !is.na(s) & !is.na(n) & n >= 2
  sqrt(sum((n[ok] - 1) * s[ok]^2) / sum(n[ok] - 1))
}
