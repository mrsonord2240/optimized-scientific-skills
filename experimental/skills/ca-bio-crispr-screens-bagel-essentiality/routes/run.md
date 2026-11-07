# BAGEL2 workflow

```bash
bash scripts/run_bagel2.sh
```

Edit the script's input paths, control and treatment names, reference sets and output directory first. It uses 30 minimum control reads, seed 42 and 1000 bootstrap iterations; it downloads missing reference files. Check available references and installed BAGEL flags beforehand. Pre/post checks guard reference overlap and NaN scores. See `references/workflow.md` and `references/errors.md`; regenerate precision-recall for the screen.

Done when fold changes, Bayes factors and precision-recall outputs pass checks and are recorded.
