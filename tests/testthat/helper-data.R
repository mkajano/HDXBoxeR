## Synthetic HDXExaminer "All results" exports for tests.
##
## make_hdx_csv() writes a csv with the same header as the shipped example,
## with deterministic uptake values so tests can check exact numbers:
##   uptake(state s, time t, peptide p, replicate r) =
##     base[p] * time_factor[t] + state_shift[s] + 0.01 * r
## Percent deuteration = 100 * uptake / nb_exch_deut(peptide).
## Undeuterated ("0s") and fully deuterated ("FD") rows are added with
## Exp Cent values (m/z) so back-exchange can be computed:
##   Exp Cent(0s) = mono_mz, Exp Cent(FD) = mono_mz + fd_uptake / charge.

hdx_header <- c("Protein State", "Deut Time", "Experiment", "Start", "End", "Sequence",
                "Adduct", "Charge", "Search RT", "Actual RT", "# Spectra", "Search IMS",
                "IMS Range", "Peak Width", "m/z Shift", "Max Inty", "Exp Cent", "Theor Cent",
                "Score", "Cent Diff", "# Deut", "Deut %", "Mark", "Confidence")

default_peptides <- data.frame(
  Start    = c(1, 5, 12, 20),
  End      = c(10, 14, 19, 30),
  Sequence = c("MKTAYIAKQR", "YIAKQRPISF", "PISFVKSH", "FSRQLEERLGL"),
  Charge   = c(2, 2, 1, 3),
  RT       = c(5.1, 6.2, 7.3, 8.4),
  base     = c(3, 4, 2, 5),
  stringsAsFactors = FALSE
)

make_hdx_csv <- function(states = c("A", "B"),
                         times = c("10.00s", "100.00s"),
                         replicates = 3,
                         peptides = default_peptides,
                         state_shift = NULL,
                         time_factor = NULL,
                         drop = NULL,
                         controls = TRUE,
                         fd_uptake = NULL,
                         path = tempfile(fileext = ".csv")) {
  if (length(replicates) == 1) replicates <- stats::setNames(rep(replicates, length(states)), states)
  if (is.null(names(replicates))) names(replicates) <- states
  if (is.null(state_shift)) state_shift <- stats::setNames((seq_along(states) - 1) * 0.5, states)
  if (is.null(time_factor)) time_factor <- stats::setNames(seq(0.3, 0.9, length.out = length(times)), times)
  nb <- HDXBoxeR::nb_exch_deut(peptides)
  mono_mz <- 400 + 50 * seq_len(nrow(peptides))
  if (is.null(fd_uptake)) fd_uptake <- 0.7 * nb

  row <- function(state, time, exp, p, uptake, cent) {
    data.frame(state, time, exp, peptides$Start[p], peptides$End[p], peptides$Sequence[p],
               "H", peptides$Charge[p], peptides$RT[p], "", 2, "", "", 0.05, 0.05, 1e5,
               cent, mono_mz[p], 0.9, "",
               if (is.na(uptake)) "" else uptake,
               if (is.na(uptake)) "" else 100 * uptake / nb[p],
               "", "High", stringsAsFactors = FALSE)
  }

  rows <- list()
  for (s in states) {
    if (controls) {
      for (p in seq_len(nrow(peptides))) {
        rows[[length(rows) + 1]] <- row(s, "0s", paste0(s, "_1"), p, NA, mono_mz[p])
        rows[[length(rows) + 1]] <- row(s, "FD", paste0(s, "_FD"), p, fd_uptake[p],
                                        mono_mz[p] + fd_uptake[p] / peptides$Charge[p])
      }
    }
    for (t in times) {
      for (r in seq_len(replicates[[s]])) {
        for (p in seq_len(nrow(peptides))) {
          u <- peptides$base[p] * time_factor[[t]] + state_shift[[s]] + 0.01 * r
          rows[[length(rows) + 1]] <- row(s, t, paste0(s, "_", t, "_", r), p, u,
                                          mono_mz[p] + u / peptides$Charge[p])
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  colnames(out) <- hdx_header
  if (!is.null(drop)) {
    ## drop: data.frame with columns state, time, replicate, Start
    for (k in seq_len(nrow(drop))) {
      hit <- out[["Protein State"]] == drop$state[k] & out[["Deut Time"]] == drop$time[k] &
        out$Experiment == paste0(drop$state[k], "_", drop$time[k], "_", drop$replicate[k]) &
        out$Start == drop$Start[k]
      out <- out[!hit, ]
    }
  }
  utils::write.csv(out, path, row.names = FALSE)
  path
}

## expected uptake for one cell of make_hdx_csv() defaults
expected_uptake <- function(s_idx, t_idx, p_idx, r, n_times = 2, n_states = 2,
                            peptides = default_peptides) {
  tf <- seq(0.3, 0.9, length.out = n_times)[t_idx]
  peptides$base[p_idx] * tf + (s_idx - 1) * 0.5 + 0.01 * r
}
