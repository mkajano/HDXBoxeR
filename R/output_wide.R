## Shared engine behind output_tp() and output_tc(): long HDXExaminer rows ->
## wide table with one column per group x replicate.
##
## layout = "tp": rows are peptide x Deut.Time, groups are protein states.
## layout = "tc": rows are peptide x Protein.State, groups are Deut.Times.
##
## Replicates: within each state and time, experiments are numbered in file
## order. All of them are used unless `replicates` caps the number. A group's
## width is its largest replicate count; a missing measurement is NA.
## A row is kept only if every group has at least `min_replicates` values.
## seq_match = FALSE: peptides must have the same Sequence in every group.
## seq_match = TRUE: peptides are matched by position (Start, End, Search.RT,
## Charge); differing sequences are joined with "." in the Sequence column.

.hdx_wide <- function(a, states, times, replicates = NULL, min_replicates = 2,
                      seq_match = FALSE, percent = FALSE, layout = c("tp", "tc")) {
  layout <- match.arg(layout)
  value_col <- if (isTRUE(percent)) "Deut.." else "X..Deut"
  value_tag <- value_col

  a <- a[a$Protein.State %in% states & a$Deut.Time %in% times, , drop = FALSE]
  a$.value <- suppressWarnings(as.numeric(a[[value_col]]))
  a <- a[!is.na(a$.value), , drop = FALSE]

  ## replicate number within state x time, in file order of Experiment
  st_time <- paste(a$Protein.State, a$Deut.Time, sep = "\r")
  a$.rep <- stats::ave(seq_len(nrow(a)), st_time, FUN = function(i) match(a$Experiment[i], unique(a$Experiment[i])))
  if (!is.null(replicates)) a <- a[a$.rep <= replicates, , drop = FALSE]

  if (layout == "tp") {
    row_vars <- c("Deut.Time", "Start", "End", "Search.RT", "Charge")
    group_var <- "Protein.State"; group_levels <- states
    group_label <- function(g) paste0(gsub(" ", "", g), "_", value_tag)
  } else {
    row_vars <- c("Protein.State", "Start", "End", "Search.RT", "Charge")
    group_var <- "Deut.Time"; group_levels <- times
    group_label <- function(g) paste0("t", g, "_", value_tag)
  }
  if (!isTRUE(seq_match)) row_vars <- c(row_vars, "Sequence")

  row_key <- do.call(paste, c(a[row_vars], sep = "\r"))
  grp <- a[[group_var]]

  ## one value per row x group x replicate; keep the first if duplicated
  cell <- paste(row_key, grp, a$.rep, sep = "\n")
  dup <- duplicated(cell)
  if (any(dup)) {
    warning(sprintf("%d duplicated measurement(s) (same peptide, state, time and replicate) - the first one is used",
                    sum(dup)), call. = FALSE)
    a <- a[!dup, , drop = FALSE]; row_key <- row_key[!dup]; grp <- grp[!dup]
  }

  ## only rows present in every group can be compared
  keys <- unique(row_key)
  width <- vapply(group_levels, function(g) {
    r <- a$.rep[grp == g]; if (length(r)) max(r) else 0L
  }, numeric(1))
  if (any(width == 0)) {
    stop(sprintf("no data for: %s", paste(group_levels[width == 0], collapse = ", ")), call. = FALSE)
  }

  ki <- match(row_key, keys)
  blocks <- lapply(seq_along(group_levels), function(j) {
    m <- matrix(NA_real_, nrow = length(keys), ncol = width[j])
    sel <- grp == group_levels[j]
    m[cbind(ki[sel], a$.rep[sel])] <- a$.value[sel]
    m
  })
  n_per_group <- sapply(blocks, function(m) rowSums(!is.na(m)))
  if (is.null(dim(n_per_group))) n_per_group <- matrix(n_per_group, nrow = 1)
  keep <- apply(n_per_group >= min_replicates, 1, all)

  ## descriptors of each kept row: first occurrence; Sequence per group
  first <- match(keys, row_key)
  desc <- a[first, c("Protein.State", "Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge")]
  if (isTRUE(seq_match)) {
    seq_by_group <- sapply(group_levels, function(g) {
      sel <- grp == g
      a$Sequence[sel][match(keys, row_key[sel])]
    })
    if (is.null(dim(seq_by_group))) seq_by_group <- matrix(seq_by_group, nrow = 1)
    desc$Sequence <- apply(seq_by_group, 1, function(s) {
      s <- s[!is.na(s)]
      if (length(unique(s)) <= 1) s[1] else paste(s, collapse = ".")
    })
  }

  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    message(sprintf("%d peptide row(s) dropped: fewer than %d replicate value(s) in at least one %s",
                    n_dropped, min_replicates, if (layout == "tp") "state" else "time point"))
  }

  value_df <- do.call(cbind, lapply(seq_along(group_levels), function(j) {
    m <- blocks[[j]][keep, , drop = FALSE]
    colnames(m) <- paste0(group_label(group_levels[j]), "_", seq_len(width[j]))
    m
  }))
  desc <- desc[keep, , drop = FALSE]

  if (layout == "tp") {
    out <- data.frame(desc[, c("Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge")],
                      value_df, check.names = TRUE)
    t_ord <- suppressWarnings(as.numeric(str_sub(out$Deut.Time, end = -2)))
    o <- order(t_ord, out$Deut.Time, out$Start, out$End, out$Charge)
  } else {
    out <- data.frame(desc[, c("Sequence", "Protein.State", "Start", "End", "Search.RT", "Charge")],
                      value_df, check.names = TRUE)
    o <- order(out$Start, out$End, out$Charge, match(out$Protein.State, states))
  }
  out <- out[o, , drop = FALSE]
  rownames(out) <- NULL

  reps <- stats::setNames(as.integer(width), group_levels)
  attr(out, "replicates") <- reps
  attr(out, "states") <- states
  attr(out, "times") <- times
  attr(out, "percent") <- isTRUE(percent)
  out
}

## message helper: one line, no stray NULL from message(cat(...))
.hdx_msg <- function(label, x) message(paste(label, paste(x, collapse = " ")))
