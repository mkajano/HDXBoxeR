## Input layer. Every function that takes `filepath` reads it through
## .hdx_raw(), so `filepath` may be a path to an HDXExaminer or DynamX csv,
## or a table returned by read_hdx() / as_hdx_long().
##
## The internal format is HDXExaminer's "All results" table, with R's
## column names: Protein.State, Deut.Time ("3.00s"; "0s" undeuterated,
## "FD" fully deuterated), Experiment (replicate), Start, End, Sequence,
## Charge, Search.RT, X..Deut (uptake, Da), Deut.. (percent), Exp.Cent
## (centroid m/z, used for back-exchange).

.hdx_required <- c("Protein.State", "Deut.Time", "Experiment", "Start", "End", "Sequence",
                   "Charge", "Search.RT", "X..Deut", "Deut..")
.proton_mass <- 1.00727646688

## the whole input table, from a path or an already-read table
.hdx_raw <- function(x) {
  if (is.data.frame(x)) {
    miss <- setdiff(.hdx_required, colnames(x))
    if (length(miss)) stop("input table missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
    if (!"Exp.Cent" %in% colnames(x)) x$Exp.Cent <- NA_real_
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  if (!is.character(x) || length(x) != 1) {
    stop("filepath must be a file path or a table from read_hdx()/as_hdx_long()", call. = FALSE)
  }
  fmt <- .hdx_detect_format(x)
  if (fmt == "hdxexaminer") return(read_hdxexaminer(x))
  read_dynamx(x, type = fmt)
}

.hdx_detect_format <- function(path) {
  if (!file.exists(path)) stop("file not found: ", path, call. = FALSE)
  first <- readLines(path, n = 1, warn = FALSE)
  first <- sub("^\ufeff", "", first)
  hdr <- trimws(strsplit(first, ",", fixed = TRUE)[[1]])
  hdr <- gsub('^"|"$', "", hdr)
  if ("Protein State" %in% hdr || "Protein.State" %in% hdr) return("hdxexaminer")
  if (all(c("Exposure", "Center") %in% hdr)) {
    return(if (all(c("File", "z") %in% hdr)) "cluster" else "state")
  }
  stop("unrecognised input format (expected an HDXExaminer 'All results' table ",
       "or a DynamX cluster/state export): ", path, call. = FALSE)
}

#' Read HDX-MS data from HDXExaminer or DynamX
#'
#' Reads an HDX-MS export into the table every HDXBoxeR function works on.
#' The result can be passed wherever a function asks for `filepath`, e.g.
#' `output_tp(read_hdx("cluster.csv", fd = "Full Deuteration control"))`.
#' Passing a file path directly also works; reading once is faster when
#' several functions use the same data, and is needed for DynamX options
#' such as `fd`.
#'
#' Supported formats (detected from the header line):
#' \describe{
#'   \item{HDXExaminer}{"All results" table export ([read_hdxexaminer()]).}
#'   \item{DynamX cluster data}{one row per peptide, state, exposure, raw
#'     file and charge ([read_dynamx()]); DynamX 2.x and 3.x. Raw files are
#'     the replicates.}
#'   \item{DynamX state data}{one averaged row per peptide, state and
#'     exposure. It holds no replicates, so it can be plotted but the
#'     replicate-based significance tests need cluster data.}
#' }
#' Any other software: build the table from a long data frame with
#' [as_hdx_long()].
#'
#' @param path path to the csv file.
#' @param format `"auto"` (default) detects the format from the header;
#'   or one of `"hdxexaminer"`, `"dynamx_cluster"`, `"dynamx_state"`.
#' @param ... passed to [read_dynamx()] (e.g. `fd`, `protein`,
#'   `combine_charges`).
#' @return data frame of class `hdx_data` in HDXExaminer's column layout.
#' @seealso [read_dynamx()], [as_hdx_long()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' d <- read_hdx(file_nm)
#' a <- output_tp(d)
#' @export
read_hdx <- function(path, format = c("auto", "hdxexaminer", "dynamx_cluster", "dynamx_state"), ...) {
  format <- match.arg(format)
  if (format == "auto") {
    fmt <- .hdx_detect_format(path)
    format <- switch(fmt, hdxexaminer = "hdxexaminer", cluster = "dynamx_cluster", state = "dynamx_state")
  }
  switch(format,
         hdxexaminer = read_hdxexaminer(path),
         dynamx_cluster = read_dynamx(path, type = "cluster", ...),
         dynamx_state = read_dynamx(path, type = "state", ...))
}

#' Read an HDXExaminer "All results" table
#'
#' Reads the csv exported by HDXExaminer (Tools > Export > All results
#' table), including exports whose data rows carry one more field than the
#' header (the "Confidence column shift").
#'
#' @param path path to the csv file.
#' @return data frame of class `hdx_data`.
#' @seealso [read_hdx()]
#' @examples
#' file_nm <- system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")
#' d <- read_hdxexaminer(file_nm)
#' @export
read_hdxexaminer <- function(path) {
  a <- utils::read.csv(file = path, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  nm <- utils::read.csv(file = path, header = TRUE, row.names = NULL, nrows = 1)
  nm1 <- colnames(nm)
  if (length(unique(nm1 == "row.names")) == 2) {
    nm1 <- nm1[-which(nm1 == "row.names")]
    nm1 <- c(nm1, "row.names")
  }
  colnames(a) <- nm1
  miss <- setdiff(.hdx_required, nm1)
  if (length(miss)) stop("input file missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  if (!"Exp.Cent" %in% colnames(a)) a$Exp.Cent <- NA_real_
  class(a) <- c("hdx_data", "data.frame")
  a
}

## theoretical maximum uptake as HDXExaminer ("Deut %") and DynamX
## (MaxUptake) count it: length - 1 - prolines after the first residue.
## Modification suffixes ("SEQ+mod") are ignored.
.hdx_theoretical_max <- function(sequence) {
  s <- sub("\\+.*$", "", sequence)
  nchar(s) - 1 - stringr::str_count(substr(s, 2, nchar(s)), "P")
}

## seconds -> HDXExaminer-style time label. Times of 1 s or more are rounded
## to whole seconds: DynamX stores minutes with 3 decimals (10 s = 0.167 min
## = 10.02 s), which would otherwise give labels like "10.02s".
.hdx_time_label <- function(sec) {
  sec <- ifelse(sec >= 1, round(sec), sec)
  sprintf("%.2fs", sec)
}

#' Read a DynamX cluster or state data export
#'
#' Converts a Waters DynamX export to the HDXBoxeR input table.
#'
#' **Cluster data** (one row per peptide, state, exposure, raw file and
#' charge; DynamX 2.x or 3.x): each raw `File` is a replicate. `Center` is
#' the centroid m/z, so the neutral mass is `z * (Center - proton mass)`.
#' Uptake of a replicate is its mass minus the mean mass of the
#' undeuterated (`Exposure` 0) measurements of the same peptide, charge and
#' state. Percent deuteration is `100 * uptake / MaxUptake`.
#'
#' **State data** (one averaged row per peptide, state and exposure):
#' `Uptake` is used directly as a single measurement per time point. It
#' contains no replicates, so functions that need replicate SDs or t-tests
#' cannot use it; a warning says so.
#'
#' Exposure is in minutes in DynamX and is converted to seconds (time
#' labels such as `"600.00s"`; times of 1 s or more are rounded to whole
#' seconds, so 0.167 min becomes `"10.00s"`). Very short "in-exchange"
#' controls (e.g. Exposure 0.001 min) are kept as a time point; leave them
#' out with the `times` argument of [output_tp()]. Note that DynamX's `MaxUptake` counts
#' exchangeable amides as length - 1 - prolines, while HDXBoxeR's
#' [nb_exch_deut()] uses length - 2 - prolines.
#'
#' @param path path to the csv file.
#' @param type `"auto"` (default), `"cluster"` or `"state"`.
#' @param fd how the fully deuterated control is identified: `NULL`
#'   (default, none), the name of a DynamX State holding the control (its
#'   measurements become the `"FD"` control of every other state), or
#'   `"last"` to use each state's longest exposure as its control (that
#'   exposure is then not kept as a time point).
#' @param protein which DynamX `Protein` to keep when the file holds more
#'   than one. Default: all (an error if there are several).
#' @param combine_charges `FALSE` (default) keeps every charge state as its
#'   own peptide, as for HDXExaminer input. `TRUE` combines charge states
#'   of a peptide into one intensity-weighted mass per raw file (as in the
#'   DynamX state export).
#' @return data frame of class `hdx_data`.
#' @seealso [read_hdx()]
#' @export
read_dynamx <- function(path, type = c("auto", "cluster", "state"), fd = NULL,
                        protein = NULL, combine_charges = FALSE) {
  type <- match.arg(type)
  if (type == "auto") type <- .hdx_detect_format(path)
  if (!type %in% c("cluster", "state")) stop("not a DynamX export: ", path, call. = FALSE)
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
  colnames(d) <- sub("^\ufeff", "", colnames(d))
  if ("Max Exchangers" %in% colnames(d)) colnames(d)[colnames(d) == "Max Exchangers"] <- "MaxUptake"
  if (!"Modification" %in% colnames(d)) d$Modification <- ""
  d$Modification[is.na(d$Modification)] <- ""

  need <- c("Protein", "Start", "End", "Sequence", "MaxUptake", "State", "Exposure", "Center",
            if (type == "cluster") c("File", "z", "RT", "Inten") else c("Uptake", "RT"))
  miss <- setdiff(need, colnames(d))
  if (length(miss)) stop("DynamX file missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  if (!is.null(protein)) d <- d[d$Protein %in% protein, , drop = FALSE]
  if (length(unique(d$Protein)) > 1) {
    stop("file holds several proteins (", paste(unique(d$Protein), collapse = ", "),
         "); choose one with protein = ", call. = FALSE)
  }
  d$Exposure <- round(as.numeric(d$Exposure), 3)
  d$Sequence <- ifelse(d$Modification != "", paste0(d$Sequence, "+", d$Modification), d$Sequence)

  out <- if (type == "cluster") .dynamx_cluster(d, combine_charges) else .dynamx_state(d)
  out <- .dynamx_fd(out, fd)
  out$Deut.Time[out$Deut.Time == .hdx_time_label(0)] <- "0s"
  class(out) <- c("hdx_data", "data.frame")
  out
}

.dynamx_cluster <- function(d, combine_charges) {
  d$z <- as.integer(d$z)
  d$mass <- d$z * (as.numeric(d$Center) - .proton_mass)
  d$Inten <- as.numeric(d$Inten)
  if (combine_charges) {
    key <- c("State", "Start", "End", "Sequence", "MaxUptake", "Exposure", "File")
    g <- interaction(d[key], drop = TRUE, lex.order = TRUE)
    first <- !duplicated(g)
    agg <- d[first, key]
    agg$mass <- as.numeric(tapply(d$mass * d$Inten, g, sum) / tapply(d$Inten, g, sum))[as.integer(g[first])]
    agg$RT <- as.numeric(tapply(as.numeric(d$RT), g, stats::median))[as.integer(g[first])]
    agg$z <- 0L
    agg$Center <- NA_real_
    d <- agg
  }
  pep <- c("State", "Start", "End", "Sequence", "z")
  pkey <- do.call(paste, c(d[pep], sep = "\r"))
  ud <- d$Exposure == 0
  m0 <- tapply(d$mass[ud], pkey[ud], mean)
  d$uptake <- d$mass - as.numeric(m0[pkey])
  ## a state without its own undeuterated run (e.g. an FD control) uses the
  ## undeuterated mass of the same peptide and charge from the other states;
  ## that mass does not depend on the state
  akey <- do.call(paste, c(d[c("Start", "End", "Sequence", "z")], sep = "\r"))
  m0_any <- tapply(d$mass[ud], akey[ud], mean)
  borrow <- is.na(d$uptake) & !ud & !is.na(m0_any[akey])
  if (any(borrow)) {
    message("undeuterated mass taken from the other states for: ",
            paste(unique(d$State[borrow]), collapse = ", "),
            " (", sum(borrow), " measurement(s))")
    d$uptake[borrow] <- d$mass[borrow] - as.numeric(m0_any[akey[borrow]])
  }
  ## one retention time per peptide and charge across all states, since
  ## output_tp() matches peptides between states on it
  rtkey <- do.call(paste, c(d[c("Start", "End", "Sequence", "z")], sep = "\r"))
  rt <- tapply(as.numeric(d$RT), rtkey, stats::median)

  out <- data.frame(
    Protein.State = d$State,
    Deut.Time = .hdx_time_label(d$Exposure * 60),
    Experiment = d$File,
    Start = as.integer(d$Start), End = as.integer(d$End), Sequence = d$Sequence,
    Charge = d$z,
    Search.RT = round(as.numeric(rt[rtkey]), 3),
    X..Deut = ifelse(ud, NA_real_, d$uptake),
    Deut.. = ifelse(ud, NA_real_, 100 * d$uptake / as.numeric(d$MaxUptake)),
    Exp.Cent = as.numeric(d$Center),
    stringsAsFactors = FALSE
  )
  dropped <- is.na(out$X..Deut) & !ud
  if (any(dropped)) {
    message(sum(dropped), " measurement(s) without an undeuterated (Exposure 0) reference were dropped")
    out <- out[!dropped, , drop = FALSE]
  }
  out
}

.dynamx_state <- function(d) {
  warning("DynamX state data holds one averaged value per peptide and time (no replicates): ",
          "use cluster data for significance tests", call. = FALSE)
  ud <- d$Exposure == 0
  rtkey <- do.call(paste, c(d[c("Start", "End", "Sequence")], sep = "\r"))
  rt <- tapply(as.numeric(d$RT), rtkey, stats::median)
  data.frame(
    Protein.State = d$State,
    Deut.Time = .hdx_time_label(d$Exposure * 60),
    Experiment = paste0(d$State, "_DynamX_state"),
    Start = as.integer(d$Start), End = as.integer(d$End), Sequence = d$Sequence,
    Charge = 1L,
    Search.RT = round(as.numeric(rt[rtkey]), 3),
    X..Deut = ifelse(ud, NA_real_, as.numeric(d$Uptake)),
    Deut.. = ifelse(ud, NA_real_, 100 * as.numeric(d$Uptake) / as.numeric(d$MaxUptake)),
    ## Center is MH+ (Da); store the neutral mass so (FD - 0s) x Charge(1) is a mass difference
    Exp.Cent = as.numeric(d$Center) - .proton_mass,
    stringsAsFactors = FALSE
  )
}

.dynamx_fd <- function(out, fd) {
  if (is.null(fd)) return(out)
  if (identical(fd, "last")) {
    last_by_state <- vapply(unique(out$Protein.State), function(s) {
      sel <- out$Protein.State == s & out$Deut.Time != .hdx_time_label(0)
      max(as.numeric(sub("s$", "", out$Deut.Time[sel])))
    }, numeric(1))
    if (length(unique(last_by_state)) > 1) {
      warning("fd = \"last\": the states have different longest exposures (",
              paste(sprintf("%s: %gs", names(last_by_state), last_by_state), collapse = ", "),
              "); each state's own longest exposure is used as its FD control", call. = FALSE)
    }
    for (s in unique(out$Protein.State)) {
      sel <- out$Protein.State == s & out$Deut.Time != .hdx_time_label(0)
      tsec <- as.numeric(sub("s$", "", out$Deut.Time[sel]))
      last <- .hdx_time_label(max(tsec))
      out$Deut.Time[sel & out$Deut.Time == last] <- "FD"
    }
    return(out)
  }
  if (!fd %in% out$Protein.State) {
    stop("fd state '", fd, "' not found; states are: ",
         paste(unique(out$Protein.State), collapse = ", "), call. = FALSE)
  }
  fd_rows <- out[out$Protein.State == fd & out$Deut.Time != .hdx_time_label(0), , drop = FALSE]
  fd_rows$Deut.Time <- "FD"
  rest <- out[out$Protein.State != fd, , drop = FALSE]
  copies <- lapply(unique(rest$Protein.State), function(s) { r <- fd_rows; r$Protein.State <- s; r })
  rbind(rest, do.call(rbind, copies))
}

#' Build HDXBoxeR input from a long table
#'
#' For data from any software: give one row per measurement (peptide,
#' charge, state, time, replicate) and name the columns that hold each
#' quantity. The result can be passed to any function that takes
#' `filepath`.
#'
#' @param data data frame with one row per measurement.
#' @param state,time,replicate,start,end,sequence,uptake names of the
#'   columns holding the protein state, labeling time, replicate
#'   (experiment/run name), peptide start and end residue, sequence and
#'   deuterium uptake in Da.
#' @param charge optional column with the charge state (default: every
#'   peptide charge 1).
#' @param percent optional column with percent deuteration; by default
#'   `100 * uptake / (length - 1 - prolines after the first residue)`, the
#'   theoretical maximum HDXExaminer and DynamX use.
#' @param rt optional column with retention time (default 0).
#' @param centroid optional column with the measured centroid m/z (needed
#'   only for [backHX_calculations()]).
#' @param time_unit unit of a numeric `time` column: `"s"` (default),
#'   `"min"` or `"h"`. Character times are used as they are.
#' @param undeuterated time value of undeuterated controls (default `0`);
#'   these rows become `"0s"`.
#' @param fd optional time value marking fully deuterated controls (e.g.
#'   `"FD"` or `Inf`); these rows become `"FD"`.
#' @return data frame of class `hdx_data`.
#' @seealso [read_hdx()]
#' @examples
#' long <- data.frame(state = rep(c("apo", "holo"), each = 6),
#'                    t_min = rep(rep(c(1, 10), each = 3), 2),
#'                    run = rep(paste0("r", 1:3), 4),
#'                    from = 5, to = 14, seq = "YIAKQRPISF",
#'                    D = c(1.1, 1.2, 1.15, 2.4, 2.5, 2.45, 0.9, 1.0, 0.95, 2.0, 2.1, 2.05))
#' d <- as_hdx_long(long, state = "state", time = "t_min", replicate = "run",
#'                  start = "from", end = "to", sequence = "seq", uptake = "D",
#'                  time_unit = "min")
#' output_tp(d)
#' @export
as_hdx_long <- function(data, state, time, replicate, start, end, sequence, uptake,
                        charge = NULL, percent = NULL, rt = NULL, centroid = NULL,
                        time_unit = c("s", "min", "h"), undeuterated = 0, fd = NULL) {
  time_unit <- match.arg(time_unit)
  cols <- c(state, time, replicate, start, end, sequence, uptake, charge, percent, rt, centroid)
  miss <- setdiff(cols, colnames(data))
  if (length(miss)) stop("column(s) not found in data: ", paste(miss, collapse = ", "), call. = FALSE)

  tm <- data[[time]]
  is_ud <- !is.na(tm) & as.character(tm) == as.character(undeuterated)
  is_fd <- if (is.null(fd)) rep(FALSE, length(tm)) else !is.na(tm) & as.character(tm) == as.character(fd)
  lab <- if (is.numeric(tm)) {
    .hdx_time_label(tm * c(s = 1, min = 60, h = 3600)[[time_unit]])
  } else as.character(tm)
  lab[is_ud] <- "0s"; lab[is_fd] <- "FD"

  out <- data.frame(
    Protein.State = as.character(data[[state]]),
    Deut.Time = lab,
    Experiment = as.character(data[[replicate]]),
    Start = as.integer(data[[start]]), End = as.integer(data[[end]]),
    Sequence = as.character(data[[sequence]]),
    Charge = if (is.null(charge)) 1L else as.integer(data[[charge]]),
    Search.RT = if (is.null(rt)) 0 else as.numeric(data[[rt]]),
    X..Deut = as.numeric(data[[uptake]]),
    stringsAsFactors = FALSE
  )
  out$Deut.. <- if (is.null(percent)) 100 * out$X..Deut / .hdx_theoretical_max(out$Sequence) else as.numeric(data[[percent]])
  out$Exp.Cent <- if (is.null(centroid)) NA_real_ else as.numeric(data[[centroid]])
  out$X..Deut[is_ud] <- NA_real_; out$Deut..[is_ud] <- NA_real_
  class(out) <- c("hdx_data", "data.frame")
  out
}
