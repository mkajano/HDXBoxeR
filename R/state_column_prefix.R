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
