# Example input: All_results_table.csv from HDXExaminer

The package ships one example HDXExaminer "All results" export, used by
the examples, tests and vignettes. It is a file, not an R object - get
its path with
`system.file("extdata", "All_results_table.csv", package = "HDXBoxeR")`.
It holds two protein states ("Unbound", "bound"), four labeling times (3
s, 60 s, 1800 s, 72000 s) plus undeuterated (`0s`) and fully deuterated
(`FD`) controls, and three replicates per state and time.

## Format

A csv file with one row per peptide, charge state, protein state,
labeling time and replicate:

- Protein State:

  protein state name

- Deut Time:

  labeling time, e.g. `"3.00s"`; `"0s"` undeuterated, `"FD"` fully
  deuterated

- Experiment:

  replicate (experiment) name

- Start:

  first residue number of the peptide

- End:

  last residue number of the peptide

- Sequence:

  amino acid sequence

- Adduct:

  adduct ion

- Charge:

  peptide charge state

- Search RT:

  search retention time

- Actual RT:

  actual retention time range

- \# Spectra:

  number of spectra

- Search IMS:

  search ion-mobility value

- IMS Range:

  ion-mobility range

- Peak Width:

  peak width

- m/z Shift:

  m/z shift

- Max Inty:

  maximum intensity

- Exp Cent:

  measured centroid (m/z)

- Theor Cent:

  theoretical centroid (m/z)

- Score:

  fit score

- Cent Diff:

  centroid difference

- \# Deut:

  deuterium uptake (Da)

- Deut %:

  percent deuteration

- Mark:

  user mark

- Confidence:

  HDXExaminer confidence level

## Source

Cook M, Delbecq SP, Schweppe TP, Guttman M, Klevit RE, Brzovic PS. The
ubiquitin ligase SspH1 from Salmonella uses a modular and dynamic E3
domain to catalyze substrate ubiquitylation. J Biol Chem.
2019;294(3):783-793.
[doi:10.1074/jbc.RA118.004247](https://doi.org/10.1074/jbc.RA118.004247)

## Details

How to prepare this export from your own HDXExaminer project:

1.  Peptides tab \> Pool \> right-click a peptide - delete peptides with
    no usable data.

2.  Tools \> Options \> Display \> Deuteration Table - display all
    columns (low-confidence results were not shown for this file).

3.  Tools \> Export \> All results table.
