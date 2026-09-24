# Protein Groups from the CLI: Percolator and Philosopher

### Protein Groups from the CLI: Percolator and Philosopher

**Goal:** Get groups and protein-level FDR out of a search you are already post-processing on the command line, and recognise a failed run instead of reporting its "0.00 % FDR".

**Approach:** Percolator does picked-protein inference itself from the `.pin` -- `-f/--picked-protein` takes the TARGET+DECOY FASTA, in-silico-digests it to build protein groups, eliminates fragment (subsumable) and duplicate proteins, then picks. The FragPipe/TPP route chains PeptideProphet -> ProteinProphet -> `philosopher filter`. Every step below is checked on the file it wrote (see the Output contract in `SKILL.md`).

```bash
# --- Percolator 3.09.0: picked-protein FDR, no extra tool after Percolator ---
# Fido (the old `--protein` / `-A`) is GONE in 3.09.0; `-f` is the only protein route.
# -P is the decoy prefix in the FASTA; -z must match the search enzyme (default trypsin).
# Add --protein-report-duplicates --protein-report-fragments (-g -c) to list group members.
DECOY_PREFIX=DECOY_   # tool-specific: DECOY_ OpenMS/Comet, rev_ Philosopher, REV__ MaxQuant
percolator -f target_decoy.fasta -P "$DECOY_PREFIX" -z trypsin \
           -l prot.target.tsv -L prot.decoy.tsv \
           -r pep.target.tsv  -B pep.decoy.tsv \
           -m psm.target.tsv  -M psm.decoy.tsv -S 1 search.pin
# OUTPUT CHECK. prot.target.tsv columns: ProteinId, ProteinGroupId, q-value,
# posterior_error_prob, peptideIds -- one row per picked representative. By default the
# indistinguishable partners are ELIMINATED, not listed: on the reference .pin all 2,339 rows
# had one accession each, and with the two report flags 28 rows listed several. Never report
# prot.target.tsv as a group list without the flags.
awk -F'\t' 'NR>1 && $3<=0.01' prot.target.tsv | wc -l    # groups at 1% protein-group FDR

# --- Philosopher 5.1.0: the FragPipe / TPP route ---
philosopher workspace --init                              # once per output folder
philosopher database --annotate target_decoy.fasta --prefix "$DECOY_PREFIX"
philosopher peptideprophet --database target_decoy.fasta --ppm --accmass \
                           --nonparam --decoy "$DECOY_PREFIX" --decoyprobs search.pep.xml
# OUTPUT CHECK: peptideprophet exits 0 even when it modelled nothing, copying every
# spectrum_query through with no probability attached. Count the results, not the exit code.
grep -c peptideprophet_result interact-search.pep.xml || { echo 'PeptideProphet modelled 0 PSMs'; exit 1; }
philosopher proteinprophet --maxppmdiff 2000000 interact-search.pep.xml
test -s interact.prot.xml || { echo 'ProteinProphet wrote no prot.xml'; exit 1; }
# --tag defaults to 'rev_', NOT the OpenMS/Comet 'DECOY_'; --razor is silently ignored
# unless --protxml supplies inference data.
philosopher filter --psm 0.01 --pep 0.01 --prot 0.01 --tag "$DECOY_PREFIX" --picked --razor \
                   --pepxml interact-search.pep.xml --protxml interact.prot.xml
philosopher report
# OUTPUT CHECK: `filter` exits 0 and prints "Converged to 0.00 % FDR with 0 PSMs" on an
# empty read, and writes header-only tables or none at all.
[ -s protein.tsv ] && [ "$(wc -l < protein.tsv)" -gt 1 ] || { echo 'filter converged on an EMPTY result'; exit 1; }
```
