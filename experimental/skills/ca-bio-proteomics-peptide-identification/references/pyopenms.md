# Database search and FDR with pyOpenMS

Read when searching or filtering inside Python with pyOpenMS 3.5 instead of an external engine.

### Database Search with pyOpenMS

**Goal:** Match tandem mass spectra in an mzML file against a protein FASTA and produce scored PSMs as idXML.

**Approach:** `SimpleSearchEngineAlgorithm` actually scores spectra (the hand-rolled `ProteaseDigestion` loop only digests, it never matches a spectrum). The FASTA must already contain target + decoy sequences concatenated for downstream FDR; decoys carry a recognizable prefix (or set `decoys` to `'true'` to let the engine generate them). Defaults are 10 ppm fragment tolerance and 1 missed cleavage, so set parameters explicitly. The search already annotates `target_decoy` on each hit.

```bash
python scripts/pyopenms_search.py sample.mzML human_target_decoy.fasta search_results.idXML \
  --precursor-ppm 10 --fragment-da 0.02 --missed-cleavages 2   # writes protein_ids FIRST, then peptide_ids
```

`scripts/pyopenms_search.py` builds the `PeptideIdentificationList` (a plain `[]` raises `TypeError` on pyOpenMS 3.5+), sets Carbamidomethyl (C) fixed and Oxidation (M) variable, and stores the idXML.

### Annotate Target/Decoy and Estimate FDR with pyOpenMS

**Goal:** Convert raw PSM scores into q-values and keep only PSMs at 1% FDR.

**Approach:** `PeptideIndexing` maps each PSM back to proteins and flags target vs decoy from the decoy prefix -- needed for idXML from other engines or after changing the FASTA; `SimpleSearchEngineAlgorithm` output above is already annotated and can go straight to `FalseDiscoveryRate`. `FalseDiscoveryRate.apply` runs the concatenated competition; `IDFilter` keeps q <= 0.01. This is the real pyOpenMS path -- not a hand-rolled decoy/target ratio of unknown provenance.

```bash
python scripts/pyopenms_fdr.py search_results.idXML human_target_decoy.fasta psms_1pct.idXML \
  --decoy-string DECOY_ --fdr 0.01     # must match the decoy prefix in the FASTA
```

`scripts/pyopenms_fdr.py` runs `PeptideIndexing` (flags target/decoy on every hit), `FalseDiscoveryRate().apply` (per-PSM q-value becomes the score), `IDFilter` at q <= 0.01, drops decoy hits and removes the spectra left with no hit (without that last step `peptide_ids.size()` still counts them: 381 instead of 311 on the synthetic run).
