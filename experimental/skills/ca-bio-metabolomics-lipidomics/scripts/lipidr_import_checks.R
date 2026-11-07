# Purpose : make lipidr's importer see LIPID MAPS 2020 sphingoid names (';O1/;O2/;O3') and fail
#           loud on any lipid it still cannot classify (Class = NA drops out of every class step).
# Inputs  : a CSV with a 'Molecule' column plus one intensity column per sample.
# Usage   : Rscript lipidr_import_checks.R in.csv [out_converted.csv]
#           or  source('lipidr_import_checks.R')  # gives to_lipidr_sphingoid(), import_with_class_check()
# Checked : lipidr 2.20.0
suppressMessages(library(lipidr))

# lipidr can't parse the ';O#' sphingoid suffix -- rewrite to the 'd/m/t' prefix it expects
# (;O1 -> m, ;O2 -> d, ;O3 -> t) before as_lipidomics_experiment()/read_skyline().
to_lipidr_sphingoid <- function(x) {
  x <- sub('(\\d+:\\d+);O1\\b', 'm\\1', x)
  x <- sub('(\\d+:\\d+);O2\\b', 'd\\1', x)
  x <- sub('(\\d+:\\d+);O3\\b', 't\\1', x)
  x
}

# raw = data.frame with Molecule + one column per sample. Converts, imports, warns on Class = NA.
import_with_class_check <- function(raw, convert = TRUE) {
  if (convert) raw$Molecule <- to_lipidr_sphingoid(raw$Molecule)   # e.g. 'Cer 18:1;O2/16:0' -> 'Cer d18:1/16:0'
  d <- as_lipidomics_experiment(raw)
  bad <- rownames(d)[is.na(rowData(d)$Class)]
  if (length(bad) > 0) {
    o_hint <- grepl(';O[0-9]', rowData(d)[bad, 'Molecule'])
    warning(sprintf('%d lipid(s) have Class = NA and would drop out of every class-based step: %s%s',
                    length(bad), paste(rowData(d)[bad, 'Molecule'], collapse = ', '),
                    if (any(o_hint)) " -- ';O#' names remain: run to_lipidr_sphingoid() before import" else ''))
  }
  d
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) < 1) stop('usage: Rscript lipidr_import_checks.R in.csv [out_converted.csv]')
  raw <- read.csv(a[1], check.names = FALSE)
  d <- import_with_class_check(raw)
  cat(sprintf('Imported %d lipids x %d samples; Class = NA: %d\n', nrow(d), ncol(d),
              sum(is.na(rowData(d)$Class))))
  if (length(a) >= 2) {
    raw$Molecule <- to_lipidr_sphingoid(raw$Molecule)
    write.csv(raw, a[2], row.names = FALSE)
  }
}
