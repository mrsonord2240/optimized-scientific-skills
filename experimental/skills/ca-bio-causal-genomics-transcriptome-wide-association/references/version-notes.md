# Source compatibility notes

These historical notes have not been reverified in this draft. Inspect installed versions before acting.

## Tool Install Notes

- **FUSION**: gusevlab.org/projects/fusion or `git clone https://github.com/gusevlab/fusion_twas`. R scripts; needs plink (PLINK 1.9) and Rscript; `install.packages(c('plink2R', 'optparse', 'glmnet', 'methods', 'RColorBrewer', 'GBJ'))` (GBJ is for the omnibus test). Pre-computed weights for GTEx v7/v8, CMC, YFS, METSIM, NTR, MESA at the same site.
- **MetaXcan / S-PrediXcan / S-MultiXcan**: `git clone https://github.com/hakyimlab/MetaXcan` (not on PyPI). Python scripts in `software/`. Compatible with Python 3.9-3.11.
- **PredictDB models**: predictdb.org. GTEx v8 elastic-net (single-tissue) and MASHR (cross-tissue posterior) databases for EUR; multi-ethnic panels emerging.
- **UTMOST**: `git clone https://github.com/Joker-Jerome/UTMOST`. Python. Needs precomputed cross-tissue weights or training pipeline.
- **FOCUS**: `pip install pyfocus "pandas<2.2" "setuptools<81"` -- checked against pyfocus 0.802 (pandas 2.1.4, numpy 1.26.4, setuptools 80.10.2), 2026-09-21. A bare `pip install pyfocus` is non-functional and the pins alone are not enough; four upstream problems, the first two fixed by the pins and the last two by the patch below (unfixed upstream as of 2026-09-21):
  - `delim_whitespace`: pandas>=2.2 removed it (used in `gwas.py`/`exprref.py`/`ldref.py`/`convert.py`), hence `pandas<2.2`.
  - `pkg_resources`: imported by pyfocus but not declared as a dependency, and gone from setuptools>=81, hence `setuptools<81`.
  - `np.warnings`: removed from every numpy that has Python 3.11/3.12 wheels (`ldref.py`, `exprref.py`).
  - `DataFrame.pivot(...)` with positional args (`finemap.py:1012`) and `DataFrame.append` (`finemap.py:188`): both removed in pandas 2.0, and both crash at "Calculating PIPs", after every earlier stage has passed.

  ```bash
  PYFOCUS_DIR=$(python -c "import pyfocus, os; print(os.path.dirname(pyfocus.__file__))")
  sed -i "1i import warnings" "$PYFOCUS_DIR/data/exprref.py"
  sed -i "s/np\.warnings\.catch_warnings/warnings.catch_warnings/;s/np\.warnings\.filterwarnings/warnings.filterwarnings/" \
      "$PYFOCUS_DIR/data/ldref.py" "$PYFOCUS_DIR/data/exprref.py"
  sed -i 's/\.pivot("model_id", "attr_name", "value")/.pivot(index="model_id", columns="attr_name", values="value")/' "$PYFOCUS_DIR/finemap.py"
  sed -i 's/^\( *\)df = df\.append(null_dict, ignore_index=True)/\1df = pd.concat([df, pd.DataFrame([null_dict])], ignore_index=True)/' "$PYFOCUS_DIR/finemap.py"
  ```

  Verified 2026-09-21 from a fresh venv with this pin + patch: `focus finemap` against a real pyfocus-schema DB (built as under "FOCUS Probabilistic Fine-Mapping") runs to completion and returns `pips_pop1 = 1` for the planted true gene and 3.4e-08 for `NULL.MODEL`. CLI `focus`. Pre-built DBs for GTEx v7/v8 panels at github.com/bogdanlab/focus. `focus import` additionally needs `mygene` + `rpy2` -- see the FOCUS section.
- **MA-FOCUS**: `git clone https://github.com/mancusolab/ma-focus && cd ma-focus && pip install .` (no PyPI release; install from source). Same CLI as single-ancestry FOCUS -- `focus finemap` -- with colon-separated per-ancestry sumstats / LD / weight DBs and paired ancestry codes in `--locations`.
- **TIGAR-V2**: `git clone https://github.com/yanglab-emory/TIGAR`. Python + R hybrid; ships with example data.
- **MOSTWAS**: `git clone https://github.com/bhattacharya-a-bt/MOSTWAS`. R package; install via `devtools::install_github`.
- **EpiXcan**: Bitbucket roussoslab/epixcan. Workflow-style pipeline.
