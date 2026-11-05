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
## seq_match = TRUE: peptides are matched by position (Start, End,
## Charge); differing sequences are joined with "." in the Sequence column.

.hdx_wide <- function(a, states, times, replicates = NULL, min_replicates = 2,
                      seq_match = FALSE, percent = FALSE, layout = c("tp", "tc"),
                      fd = NULL) {
  layout <- match.arg(layout)
  fd_mode <- identical(percent, "FD")
  per_res <- identical(percent, "per_residue")
  if (!(isTRUE(percent) || identical(percent, FALSE) || fd_mode || per_res)) {
    stop("percent must be FALSE, TRUE, \"FD\" or \"per_residue\"", call. = FALSE)
  }
  value_col <- if (isTRUE(percent)) "Deut.." else "X..Deut"
  value_tag <- if (fd_mode || per_res) "Deut.." else value_col

  a <- a[a$Protein.State %in% states & a$Deut.Time %in% times, , drop = FALSE]
  a$.value <- suppressWarnings(as.numeric(a[[value_col]]))
  a <- a[!is.na(a$.value), , drop = FALSE]
  if (per_res) {
    ## uptake per exchangeable amide (Da / nb_exch_deut): removes peptide length
    n_exch <- nb_exch_deut(a)
    a$.value <- ifelse(n_exch > 0, a$.value / n_exch, NA_real_)
    a <- a[!is.na(a$.value), , drop = FALSE]
  }
  if (fd_mode) a$.fd <- .hdx_fd_lookup(a, fd, seq_match)

  ## replicate number within state x time, in file order of Experiment
  st_time <- paste(a$Protein.State, a$Deut.Time, sep = "\r")
  a$.rep <- stats::ave(seq_len(nrow(a)), st_time, FUN = function(i) match(a$Experiment[i], unique(a$Experiment[i])))
  if (!is.null(replicates)) a <- a[a$.rep <= replicates, , drop = FALSE]

  if (layout == "tp") {
    row_vars <- c("Deut.Time", "Start", "End", "Charge")
    group_var <- "Protein.State"; group_levels <- states
    group_label <- function(g) paste0(gsub(" ", "", g), "_", value_tag)
  } else {
    row_vars <- c("Protein.State", "Start", "End", "Charge")
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

  ## percent of the fully deuterated control: rows are kept or dropped on the
  ## uptake values above; a peptide without a usable FD value keeps its row
  ## with NA values
  if (fd_mode) {
    blocks <- lapply(seq_along(group_levels), function(j) {
      f <- rep(NA_real_, length(keys))
      sel <- grp == group_levels[j]
      f[ki[sel]] <- a$.fd[sel]
      100 * blocks[[j]] / f
    })
    empty <- do.call(cbind, lapply(blocks, function(m) rowSums(!is.na(m)) == 0))
    no_fd <- sum(keep & rowSums(empty) > 0)
    if (no_fd > 0) {
      message(sprintf("%d peptide row(s) have no usable fully deuterated (FD) value in at least one %s: left empty (NA)",
                      no_fd, if (layout == "tp") "state" else "time point"))
    }
  }

  ## descriptors of each kept row: first occurrence; Sequence per group
  first <- match(keys, row_key)
  desc <- a[first, c("Protein.State", "Deut.Time", "Start", "End", "Sequence", "Search.RT", "Charge")]
  ## one retention time per peptide (median over all its measurements, all
  ## states): RT is a descriptor, not part of peptide identity, and a shared
  ## value keeps later merges of tables from different states consistent
  rt_vars <- c("Start", "End", "Charge", if (!isTRUE(seq_match)) "Sequence")
  rt_key <- do.call(paste, c(a[rt_vars], sep = "\r"))
  rt_med <- tapply(suppressWarnings(as.numeric(a$Search.RT)), rt_key, stats::median, na.rm = TRUE)
  desc$Search.RT <- as.numeric(rt_med[rt_key[first]])
  if (isTRUE(seq_match)) {
    seq_by_group <- sapply(group_levels, function(g) {
      sel <- grp == g
      a$Sequence[sel][match(keys, row_key[sel])]
    })
    if (is.null(dim(seq_by_group))) seq_by_group <- matrix(seq_by_group, nrow = 1)
    desc$Sequence <- apply(seq_by_group, 1, function(s) {
      s <- s[!is.na(s)]
      if (length(unique(s)) <= 1) s[1] else paste(unique(s), collapse = ".")
    })
  }

  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    absent <- !keep & apply(n_per_group == 0, 1, any)
    unit <- if (layout == "tp") "state" else "time point"
    message(sprintf("%d peptide row(s) dropped: %d not measured in every %s, %d with fewer than %d replicate value(s)",
                    n_dropped, sum(absent), unit, n_dropped - sum(absent), min_replicates))
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
    t_ord <- .hdx_time_seconds(out$Deut.Time)
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
  attr(out, "percent") <- if (fd_mode) "FD" else if (per_res) "per_residue" else isTRUE(percent)
  out
}

## Mean fully deuterated (FD) uptake (Da) per state and peptide, from the
## "FD" rows of the raw input. Charge states are pooled: uptake in Da does
## not depend on charge. FD uptake <= 0 cannot be used and is NA.
.hdx_fd_table <- function(raw) {
  fd <- raw[raw$Deut.Time == "FD", , drop = FALSE]
  fd$X..Deut <- suppressWarnings(as.numeric(fd$X..Deut))
  fd <- fd[!is.na(fd$X..Deut), , drop = FALSE]
  if (nrow(fd) == 0) {
    return(data.frame(Protein.State = character(0), Start = integer(0), End = integer(0),
                      Sequence = character(0), fd = numeric(0)))
  }
  out <- stats::aggregate(X..Deut ~ Protein.State + Start + End + Sequence, data = fd, FUN = mean)
  names(out)[names(out) == "X..Deut"] <- "fd"
  out$fd[out$fd <= 0] <- NA_real_
  out
}

## FD value for every measurement row of `a`; with position matching
## (seq_match = TRUE) the sequence is ignored
.hdx_fd_lookup <- function(a, fd, seq_match) {
  if (is.null(fd) || nrow(fd) == 0) return(rep(NA_real_, nrow(a)))
  if (isTRUE(seq_match)) {
    fd <- stats::aggregate(fd ~ Protein.State + Start + End, data = fd, FUN = mean, na.action = stats::na.pass)
    k_fd <- paste(fd$Protein.State, fd$Start, fd$End, sep = "\r")
    k_a <- paste(a$Protein.State, a$Start, a$End, sep = "\r")
  } else {
    k_fd <- paste(fd$Protein.State, fd$Start, fd$End, fd$Sequence, sep = "\r")
    k_a <- paste(a$Protein.State, a$Start, a$End, a$Sequence, sep = "\r")
  }
  fd$fd[match(k_a, k_fd)]
}

## message helper: one line, no stray NULL from message(cat(...))
.hdx_msg <- function(label, x) message(paste(label, paste(x, collapse = " ")))

## time label to seconds, for sorting: "1800.00s", "30 s", "1min", "2 h";
## NA for labels such as "FD"
.hdx_time_seconds <- function(x) {
  x <- trimws(as.character(x))
  num <- suppressWarnings(as.numeric(sub("^([0-9.eE+-]+).*$", "\\1", x)))
  unit <- tolower(trimws(sub("^[0-9.eE+-]+", "", x)))
  mult <- ifelse(unit %in% c("", "s", "sec"), 1,
          ifelse(unit %in% c("m", "min", "mins"), 60,
          ifelse(unit %in% c("h", "hr", "hrs", "hour", "hours"), 3600, NA_real_)))
  num * mult
}
