# ePRIDICT: Chromatin-Context Prediction

Moved from SKILL.md. Read when a PRIDICT2 prediction needs to be checked against the target's
chromatin context before library synthesis, or when a pegRNA's observed efficiency is far below
its sequence-only prediction.

Mathis N et al. 2024 *Nat Biotechnol* (doi:10.1038/s41587-024-02268-2) and Mathis N et al. 2025
*Nature Protocols*: an XGBoost model trained in K562 that predicts prime editing efficiency from
**chromatin context** rather than sequence. It is designed to **complement PRIDICT2, not replace
it**: pick the pegRNA with PRIDICT2 first, then check the resulting target's endogenous
targetability with ePRIDICT. The pairing helps most in regions of lower chromatin accessibility,
which sequence-only prediction cannot see.

## Install

```bash
# Public, unauthenticated. The tool runs on Linux/macOS only -- pybigwig has no Windows wheel,
# so on Windows use WSL. Verified on epridict git HEAD ddcdba360c969469a536343dd89e29179961e367.
git clone https://github.com/Schwank-Lab/epridict.git
cd epridict

# The repo ships epridict_env.yml. It pins old versions (pandas 1.5.2, xgboost 1.7.6) and lists
# `defaults`/`anaconda` channels, which may require a ToS acceptance step; creating a dedicated
# environment from conda-forge works and is verified to load the shipped model:
conda create -n epridict -c conda-forge python=3.10 numpy pandas scipy tqdm joblib xgboost pybigwig
conda activate epridict

# Then download the ENCODE bigWigs the model reads. The `light` model needs 6 datasets (~5.3 GB);
# the `full` model needs 455 (~624 GB) -- use light unless you have the disk.
./epridict_download_encode.sh light
```

**The repo's own download script does not verify what it downloaded.** It runs
`curl -L -s -o '<file>' '<url>'` for each accession, checks only curl's return code, and never
compares the resulting size against the ENCODE manifest (`epridict_download_encode.sh`). With no
`--fail` and no `--retry`, an HTTP error body or a connection dropped mid-transfer is written to
`bigwig/<accession>.bigWig` and reported as a success. A truncated bigWig does not raise on its
own -- it silently shifts the chromatin features that drive the score. Verify before trusting a
prediction:

```bash
# Each of the 6 light-model bigWigs must actually parse and carry signal. Compare sizes against
# https://www.encodeproject.org/files/<accession>/ (the manifest 'file_size' field) and then open
# each file, which catches a truncated body the size check alone would miss:
python - <<'PY'
import os, pyBigWig
REPO = os.getcwd()
accs = sorted({l.split('_')[1] for l in open(f'{REPO}/misc/ePRIDICT_slim_model_column_names.txt')})
for acc in accs:
    p = f'{REPO}/bigwig/{acc}.bigWig'
    bw = pyBigWig.open(p)
    vals = bw.values('chr3', 44843404, 44843604)
    n = sum(1 for v in vals if v is not None)
    print(f'{acc}: {os.path.getsize(p)/1e9:.3f} GB, {len(bw.chroms())} chroms, {n}/200 values at chr3:44843504')
    bw.close()
    assert n > 0, f'{acc} parses but carries no signal at this locus'
PY
# Verified sizes: ENCFF139KZL 639391582 | ENCFF601JGK 1017107689 | ENCFF834SEY 507390951
#                 ENCFF954LGE 2006979034 | ENCFF959YJV 397853542 | ENCFF972GVB (not in manifest)
```

The script resolves `./input`, `./predictions`, `./bigwig` and `misc/` **relative to the current
directory**, so run every command below from the repository root.

## Run a prediction

```bash
# Single locus. Input is a genomic position in hg38 -- NOT a pegRNA. Locate the position first
# (PRIDICT2 reports Editing_Position for the intended edit); ePRIDICT scores the chromatin context
# around that site. Y chromosome is not supported.
python epridict_prediction.py single --chromosome chr3 --position_hg38 44843504

# Full model instead of the default light one (requires the 455-dataset download):
# python epridict_prediction.py single --chromosome chr3 --position_hg38 44843504 --use_full_model
```

Real output (light model, PRIDICT2's own worked-example locus):

```
ePRIDICT score (light model):  42.53
Percentile in context of 27625 sampled genomic locations in K562: 62.02%

Output stored in .../epridict/./predictions as chr3_44843504.csv
```

```bash
# Batch mode: a CSV with exactly two columns, `chromosome` and `position_hg38`, placed in ./input
# (input-fname is resolved there, not in the current directory). Output defaults to
# predictions/<input_filename>_output.csv; override with --output-fname.
python epridict_prediction.py batch my_variants.csv
```

The per-locus CSV holds the model's 24 input features and two results; the batch CSV holds one row
per input position:

| Column | Meaning |
|--------|---------|
| `ePRIDICT_prediction_light` (per-locus) / `ePRIDICT score (light model)` (batch) | Predicted PE efficiency, 0-100 scale |
| `percentile_light_ePRIDICT_genomewide_K562` / `ePRIDICT percentile (light model)` | Percentile against `misc/genome_wide_ePRIDICT_slim_predictions_K562.csv` |
| `Error notes` (batch) | Populated when a locus could not be scored |
| 24 `*_average_value_{100,1000,2000,5000}` columns | Chromatin signal from the 6 ENCODE tracks (DNase-seq, HDAC2, H3K4me1, H3K9me3, H3K27me3, H3K4me2) at 4 bin sizes |

## Interpreting the score

- The score is **cell-line specific (K562)**. Prediction performance may differ in other lines;
  treat it as chromatin-context evidence, not a cross-line guarantee.
- A **low percentile** at a locus whose PRIDICT2 score is high is the signature of the
  chromatin-context gap in Failure Modes ("Low pegRNA efficiency despite high PRIDICT
  prediction"): the sequence is favorable but the locus is poorly accessible. Pilot those
  pegRNAs before committing library budget.
- Report the two numbers together. ePRIDICT never overrides a PRIDICT2 design; it re-ranks it.