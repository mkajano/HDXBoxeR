## internal helper - not exported. Shared by every function that needs to
## match a user-supplied Protein.State string against output_tp()'s
## wide-format column names (compare_states_tp(), derive_significant_regions(),
## plot_fissure_vent(), interactive_explorer(), residue_regions(),
## plot_hdx_structure()).
##
## output_tp()'s uptake_seq_matchF/T (see output_tp_def.R) build each raw
## replicate/average/SD column name as "<state>_<suffix>_<n>", where
## <state> is the user's Protein.State string with only spaces removed
## (gsub(' ', '', state)). That partially-sanitized string is then handed
## to data.frame()'s default check.names = TRUE, which runs make.names() on
## every column name and silently replaces any OTHER syntactically-invalid
## character (e.g. "+", "(", ")") with ".". A state like "tubulin+b5"
## therefore surfaces in the data as a "tubulin.b5_..." column, not
## "tubulin+b5_...", even though nothing downstream ever explicitly asked
## for that substitution.
##
## Every consumer that recovers "which columns belong to state X" has to
## compare against that same fully-mangled string, not the raw state name -
## comparing against the raw name either silently matches nothing (an exact
## `==` comparison, e.g. compare_states_tp() before this helper existed -
## every delta/p-value for that state came back NA with no warning) or
## errors as "no matching column found" (a startsWith()-based comparison,
## e.g. plot_fissure_vent() before this helper existed). Reproducing
## make.names(gsub(' ', '', state)) exactly is the only reliable way to
## recover the mapping, since it is exactly what R itself already applied.
.state_column_prefix <- function(state) {
  make.names(gsub(" ", "", state, fixed = TRUE))
}

## the one average/SD column of a state, e.g. "av_<state>_X..Deut",
## "sd_<state>_X..Deut", "sd_ci_<state>_Deut..": exact match on the state
## part after removing the known prefix and value suffix
.hdx_summary_col <- function(state, cols, prefixes = c("sd_ci_", "sd_", "av_", "pv_")) {
  prefix <- .state_column_prefix(state)
  root <- cols
  for (p in prefixes) root <- sub(paste0("^", p), "", root)
  is_value <- grepl("_(X\\.\\.Deut|Deut\\.\\.)$", cols)
  root <- sub("_(X\\.\\.Deut|Deut\\.\\.)$", "", root)
  hits <- cols[is_value & root == prefix]
  ## no exact match: accept an unambiguous abbreviation ("Fresh" -> "Fresh_pH7")
  if (length(hits) == 0) hits <- cols[is_value & startsWith(root, prefix)]
  if (length(hits) != 1) {
    stop(sprintf("State '%s' must match exactly one column among: %s",
                 state, paste(cols[is_value], collapse = ", ")), call. = FALSE)
  }
  hits
}

## raw replicate columns of one state in an output_tp() table: exact match
## on "<state>_X..Deut_<k>" / "<state>_Deut.._<k>" (not a prefix match, so
## "WT" does not also pick up "WT_mut"). Any number of replicates is
## accepted; states may differ.
.hdx_state_cols <- function(state, cols) {
  prefix <- .state_column_prefix(state)
  root <- sub("_(X\\.\\.Deut|Deut\\.\\.)_[0-9]+$", "", cols)
  is_value <- grepl("_(X\\.\\.Deut|Deut\\.\\.)_[0-9]+$", cols)
  hits <- cols[is_value & root == prefix]
  if (length(hits) == 0) {
    ## no exact match: accept an unambiguous abbreviation
    cand <- unique(root[is_value & startsWith(root, prefix)])
    if (length(cand) == 1) hits <- cols[is_value & root == cand]
  }
  if (length(hits) == 0) {
    stop(sprintf("No replicate columns for state '%s' among: %s",
                 state, paste(cols[is_value], collapse = ", ")), call. = FALSE)
  }
  hits
}
