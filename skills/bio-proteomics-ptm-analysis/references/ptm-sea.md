# PTM-SEA With ssGSEA2.0

**Goal:** Score PTMsigDB kinase, perturbation and pathway signatures from one signed site-level statistic.

**Approach:** PTM-SEA is ssGSEA run on site identifiers. Each row is ONE localized site, named by its +/-7 flanking sequence plus `-p` (`AALRQLRSPRRAQAP-p`). PTMsigDB scores each signature site with its reported direction (`;u` / `;d`, which never appear in your ids). Feed a SIGNED per-site statistic (log2FC or the moderated t from the protein-adjusted model), not intensities. Get ssGSEA2.0 from `github.com/broadinstitute/ssGSEA2.0` (checked on `cae7bed`); the human flanking database is `db/ptmsigdb/v1.9.1/ptm.sig.db.all.flanking.human.v1.9.1.gmt`, with mouse and rat files beside it.

**Input boundary:** For regulatory kinase, perturbation, or pathway claims, derive `stat` only from the checked runner's published paired-global label-free or TMT result. Do not build regulatory PTM-SEA input from `proxy_adjusted_sites.csv` or `proxy_ptm_model.csv`: no-global outputs are proxy-adjusted candidates, not protein-adjusted regulation evidence.

```bash
python scripts/ptmsea.py write-gct --sites sites.tsv --out sites.gct
```

`sites.tsv` has one row per LOCALIZED site with `seq_window` (MaxQuant "Sequence window": 31 aa, 15 each side) and `stat` (the signed statistic). The script keeps the central 15 (+/-7) plus `-p` as the id; ids must be unique, so duplicates are averaged (a choice: report it).

## Checked execution for unattended runs

Do not treat the appearance of GCT files as successful completion: the ssGSEA R
process must exit cleanly, leave no process in the runner-owned scope, and produce
a parseable NES/FDR table. The following invokes the raw ssGSEA command and the
reader as one checked child command. `run_checked.py` keeps the stage directory and
receipt on every outcome, but publishes the one final CSV only after clean exit and
CSV validation.

```bash
PYTHON=python
# Caller-supplied path to the supported, version-pinned R launcher for this run.
# Audit example only: F:/OpenScience/audit-envs/mass-spec-proteomics-analyst/r.sh
# A shell wrapper must itself be callable in the current environment; on native Windows,
# supply the validated Rscript.exe/platform launcher with its required R library configured.
RSCRIPT=/absolute/path/to/validated/R-launcher
SSG=ssGSEA2.0
DB=$SSG/db/ptmsigdb/v1.9.1/ptm.sig.db.all.flanking.human.v1.9.1.gmt
RUN_ID=$($PYTHON -c 'import uuid; print(uuid.uuid4().hex)')

# FINAL's parent must permit a temporary sibling: publication copies verified bytes,
# fsyncs them, then exposes the complete file by atomic no-clobber link.
OUT=ptmsea
STAGE="$OUT/.checked-$RUN_ID"
RECEIPT="$OUT/receipt-$RUN_ID.json"
FINAL="$OUT/ptmsea_scores-$RUN_ID.csv"  # must not already exist (including as a dangling symlink)

$PYTHON scripts/run_checked.py \
  --timeout 1800 --grace 0.15 \
  --receipt "$RECEIPT" --stage-dir "$STAGE" \
  --publish-source result.csv --publish-dest "$FINAL" \
  --csv result.csv --require-columns NES,FDR --min-rows 1 \
  -- "$PYTHON" -c '
import subprocess, sys
stage, rscript, ssg, db, sites, reader = sys.argv[1:]
subprocess.run([
    rscript, f"{ssg}/ssgsea-cli.R", "-i", sites, "-o", f"{stage}/run",
    "-d", db, "-z", ssg, "-n", "rank", "-w", "0.75", "-c", "z.score",
    "-t", "area.under.RES", "-s", "NES", "-p", "1000", "-m", "10",
    "-x", "TRUE", "-e", "FALSE", "-l", "FALSE",
], check=True)
subprocess.run([
    sys.executable, reader, "read", "--prefix", f"{stage}/run",
    "--out", f"{stage}/result.csv",
], check=True)
' "$STAGE" "$RSCRIPT" "$SSG" "$DB" "$PWD/sites.gct" scripts/ptmsea.py
```

`run_checked.py` owns a new process session on POSIX and a Windows Job Object on
Windows; it cleans only that scope on timeout, nonzero worker exit, or a lingering child. A nonzero R
exit, timeout, lingering owned process, missing/malformed `result.csv`, or existing
`FINAL` leaves staged evidence and a non-success receipt without publishing a result.
The runner intentionally publishes one file only; keep the four intermediate GCTs in
the private stage and do not infer success from them. The run still writes
`-pvalues.gct` and `-combined.gct` (signature size and overlap) before the reader
creates the published NES/FDR CSV.

Checked on ssGSEA2.0 `cae7bed` with PTMsigDB v1.9.1, running the three steps above (`write-gct`, the ssgsea-cli call, `read`; originally as inline blocks, re-run 2026-09-21 as `scripts/ptmsea.py`) on a planted input: the 588 human CDK1 substrate sites of `KINASE-PSP_CDK1` given a +2 shift among 2,582 sites in all (real PTMsigDB flanking sequences, 2,000 of them from other signatures). `KINASE-PSP_CDK1` ranked 1 of 100 scored signatures (NES 44.0, FDR 0.0056, the floor at 1,000 permutations, which at least four other signatures also reached); CDK2, CDK5, CDK6 and ERK2/MAPK1 scored too because they share substrates. With the same statistics shuffled across sites, CDK1 fell to rank 37 (NES 0.47, FDR 0.95) and no signature had FDR below 0.05 (SERUM and PKCZ sat at exactly 0.05).

Only 100 of the 495 signatures had the `-m 10` overlap with those sites, so read the overlap column (`-combined.gct`) before calling a signature absent, and raise `-p` for finer p-values. PTMsigDB also holds non-phospho ids (`-ac`, `-m2`, `-m3`); `write-gct` builds `-p` ids only. Signatures that share substrates co-score: name the substrate set, not an independent kinase (see "Over-reading kinase-activity output" in SKILL.md).
