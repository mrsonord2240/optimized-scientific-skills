## Shiba for Low-Coverage / Few-Replicate Designs

**Goal:** Detect differential splicing with explicit junction-imbalance correction — addresses a known false-positive source for rMATS-style methods.

**Approach:** Shiba is a pipeline configured via YAML. Install via bioconda, write an experiment table and a config file, then run `shiba.py` (`--mame` = splicing analysis only).

```bash
conda install -c bioconda shiba

# exp.tsv (tab-separated): sample  bam  group  technology (short|long); groups named in config.yaml
# config.yaml (keys used in the checked run):
#   workdir: shiba_out          gtf: annotation.gtf        experiment_table: exp.tsv
#   unannotated: False          minimum_anchor_length: 6   minimum_intron_length: 50
#   maximum_intron_length: 500000   strand: XS            only_psi: False   only_psi_group: False
#   fdr: 0.05                   delta_psi: 0.1             reference_group: Ref   alternative_group: Alt
#   minimum_reads: 10           individual_psi: True       ttest: False      excel: False
shiba.py -p 8 --mame config.yaml
# results: shiba_out/results/splicing/PSI_<SE|FIVE|THREE|MXE|RI|AFE|ALE|MSE>.txt; columns dPSI, q, "Diff events" (Yes/No)
```

- `strand: XS` needs XS-tagged BAMs; regtools 1.0 rejects the numeric `strand: 0/1/2` values Shiba passes through.
- The GTF must contain every event type with reads (Shiba 0.8.2 raised a KeyError on a GTF with SE events only).
- dPSI = alternative - reference group. `snakeshiba.smk` lives in the package's share directory and its config needs a container key, so `shiba.py` is the entry point that ran.
- Shiba (Kubota 2025 *NAR*) is reported to outperform rMATS at n=2 vs n=2 by correcting differential mappability between inclusion and skipping junctions. Not reproduced here: on a simulation without planted junction imbalance it found 22/24 strong events at 4v4 with 2 false calls, no better than rMATS or leafcutter, and made 2 calls on a null 2v2 comparison (separate simulation). Community calibration is still emerging. Full config schema: https://sika-zheng-lab.github.io/Shiba/
