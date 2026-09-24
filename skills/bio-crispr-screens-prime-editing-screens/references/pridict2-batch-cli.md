# PRIDICT2 Batch CLI and pegRNA Library Filtering

Moved from SKILL.md. Read when running PRIDICT2 (single or batch), parsing its output, or filtering a pegRNA library.

## PRIDICT and PRIDICT2 pegRNA Efficiency Prediction

**Mathis N et al 2023 *Nat Biotechnol* 41:1151 (PRIDICT v1) / 2025 *Nat Biotechnol* 43(5):712 (PRIDICT2; published online June 2024)** developed deep-learning predictors of per-pegRNA editing efficiency. PRIDICT2 is the current state of the art.

```bash
# PRIDICT2 is invoked via CLI: pridict2_pegRNA_design.py
# Single sequence input (parens in --sequence = intended edit):
python pridict2_pegRNA_design.py single \
    --sequence-name BRCA1_c5135 \
    --sequence "AGCAGCCT(C/T)CTGAATGCCC...60nt_context" \
    --output-dir predictions/ \
    --use_5folds                                              # 5-fold ensemble averaging

# Batch input from CSV:
# PRIDICT2's --input-dir defaults to ./input, NOT the current directory -- the CSV must live there
# (or pass --input-dir explicitly). --output-dir must also already exist AND hold no .csv file at all
# before running with --summarize (it lists existing .csv files there first: a missing directory
# raises FileNotFoundError, and a directory containing any .csv raises "Output directory is not
# empty" -- so a SECOND --summarize run into the same --output-dir always fails; use a fresh dir).
mkdir -p input predictions
mv variants_to_design.csv input/
# CSV columns: sequence_name, editseq (NOT "sequence" -- see SKILL.md Failure Modes)
python pridict2_pegRNA_design.py batch \
    --input-fname variants_to_design.csv \
    --output-dir predictions/ \
    --cores 3 \
    --summarize K562                                          # takes a cell-line value ('K562' or 'HEK'); a bare flag crashes argparse
# --cores defaults to 3 and the tool documents 3 as the MAXIMUM ("Maximum 3 cores to prevent memory
# issues", pridict2_pegRNA_design.py line 1111). The value is passed straight through unvalidated, so
# a higher number is accepted silently and then risks the memory exhaustion the cap exists to avoid.

# Output: per-pegRNA predictions in predictions/<sequence_name>_pegRNA_Pridict_full.csv
# Real columns (verified against actual output, not the names above earlier drafts guessed):
#   PRIDICT2_0_editing_Score_deep_K562, PRIDICT2_0_editing_Score_deep_HEK (0-100 scale, not 0-1),
#   PBSlength, RTlength, PBSrevcomp, RTrevcomp, Spacer-Sequence, pegRNA (full assembled sequence),
#   Target-Strand, Editing_Position, Correction_Type, Correction_Length, among ~50 total columns.
```

**Loading PRIDICT2 results in Python:**

```python
import pandas as pd
from pathlib import Path

def load_pridict2_predictions(prediction_dir):
    '''Load PRIDICT2 batch outputs from prediction_dir/'''
    summary = pd.read_csv(Path(prediction_dir) / '<timestamp>_summary_K562_batch_summary.csv')
    # Real columns (verified against real PRIDICT2 2026-09 output): sequence_name,
    # PRIDICT2_0_editing_Score_deep_K562, PRIDICT2_0_editing_Score_deep_HEK, PBSrevcomp,
    # RTrevcomp, pegRNA, Target-Strand, Editing_Position, among others -- NOT the
    # "PBS"/"RTT"/"predicted_efficiency" names some earlier drafts assumed.
    return summary
```

**Key determinants of PE efficiency (Mathis 2025 PRIDICT2):**

| Feature | Effect on efficiency |
|---------|----------------------|
| PBS GC content | 40-55% optimal; high GC slows annealing |
| PBS length | 11-13 nt optimal; longer for high-GC PBS |
| RTT length | 10-20 nt typical; trade-off between coverage and processivity |
| Edit position in RTT | Closest to PBS = highest efficiency |
| Chromatin context | Dominant locus effect; H3K9me3 heterochromatin ~0.8% vs ~2.2% elsewhere |
| Cell line / Cas9 expression | Variable; piloting required |
| Cell cycle phase | S/G2 = higher efficiency |

**Critical insight from Mathis 2025:** Chromatin context is a major locus-level determinant that sequence-only predictors miss, which is why ePRIDICT is designed to be combined with PRIDICT2.0 rather than replace it -- the pairing helps most in regions of lower chromatin accessibility. Run it on the positions PRIDICT2 selected (`references/epridict-chromatin.md`); for genome-scale screens, validate predictions empirically at representative loci.

**PRIDICT2's reported `Spacer-Sequence` forces a synthetic 5'-G** for U6-promoter transcription (source: `pridict2_pegRNA_design.py`, `protospacerseq = 'G' + original_seq[...]`) -- confirmed against real output, where the reported spacer matched true genomic sequence at 19/20 positions, the sole mismatch being the 5'-most base. If your spacer's true genomic first base isn't G, the ordered oligo will still differ from genomic sequence there; this is expected, not a bug.

## Run PRIDICT2 on a Custom pegRNA Library

**Goal:** Predict editing efficiency for thousands of pegRNAs before library synthesis.

**Approach:** Build a CSV with one row per intended edit (sequence + edit notation), run PRIDICT2 in batch mode, parse the per-pegRNA efficiency summary, and filter to candidates above the chosen efficiency threshold.

```bash
# Step 1: prepare batch input CSV -- required header is "editseq", NOT "sequence"
# (see Failure Modes: a wrong header exits 0 with an empty summary file)
# PRIDICT2's --input-dir defaults to ./input, NOT the current directory -- write the CSV there
# (or pass --input-dir explicitly), and pre-create --output-dir (see Step 2).
mkdir -p input predictions
cat > input/variants.csv <<EOF
sequence_name,editseq
BRCA1_R71X,AGCAGCCT(C/T)CTGAATGCCC...
MLH1_c677,GAGCTGAGC(A/G)GAGGCTCTTGAAGC...
EOF

# Step 2: run PRIDICT2 batch (--summarize takes a cell-line value, not a bare flag;
# --output-dir must already exist and be free of .csv files -- --summarize lists .csv files there
# before the run starts, so a missing directory raises FileNotFoundError and a directory already
# holding a previous run's .csv raises "Output directory is not empty" rather than a no-op)
python pridict2_pegRNA_design.py batch \
    --input-fname variants.csv \
    --output-dir predictions/ \
    --cores 3 \
    --summarize K562
```

**If the summary file is empty (just `""`) the run still exited 0** -- see Failure Modes, "Batch run exits 0 with an empty summary file", for the three causes.

```bash
# Step 3: parse and filter. The score column is PRIDICT2_0_editing_Score_deep_K562 (0-100 scale, not
# "predicted_editing_efficiency"); 50 is a project-chosen library-inclusion cutoff; keeps top 3 per intended edit.
python scripts/filter_pridict2_summary.py predictions/<timestamp>_summary_K562_batch_summary.csv     --threshold 50 --top 3 --out peg_library_filtered.csv
```
