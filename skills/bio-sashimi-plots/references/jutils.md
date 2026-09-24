## Jutils for Tool-Agnostic Output

**Goal:** Visualize differential splicing output uniformly across rMATS, leafcutter, MntJULiP, and MAJIQ.

**Approach:** Convert tool output to Jutils' standard TSV, then plot. Run from the Jutils clone (`python3 jutils.py ...`). `scripts/jutils_pipeline.sh` below covers the rMATS path end to end.

```bash
scripts/jutils_pipeline.sh /path/to/Jutils rmats_output/ meta.tsv bam_list.tsv annotation.gtf chr1:1000-2000 jutils_out/
```

`meta.tsv` is `sample<TAB>condition` and `bam_list.tsv` is `sample<TAB>bam<TAB>condition`. The script runs, in order: `convert-results` (writes `rmats_JC_results.tsv` and `rmats_JCEC_results.tsv`), `heatmap` (`--q-value 0.05`, env `Q`; **needs >= 2 events passing the cutoffs**, otherwise it prints "Skipping"; writes `clustermap*.pdf`), `sashimi` for the coordinate, and `venn-diagram`. `--tsv-file-list` is a FILE with one `path<TAB>label` line per TSV, not a comma-separated list; the script writes it. Outputs go to `hm/`, `sh/`, `vn/` under the output directory.

**leafcutter / MntJULiP / MAJIQ converters** (2026-09-21: run end to end on Jutils' own shipped test data, `data/` in the Jutils clone — a real published mouse hippocampus epileptic-vs-control RNA-seq set, not synthetic):

```bash
python3 jutils.py convert-results --leafcutter-dir data/leafcutter --mntjulip-dir data/mntjulip \
    --majiq-dir data/majiq --rmats-dir data/rmats --out-dir jutils_out/
```

- `--leafcutter-dir` needs `leafcutter_ds_cluster_significance.txt`, `leafcutter_ds_effect_sizes.txt`, `results_perind.counts.gz`, `results_perind_numers.counts.gz` (leafcutter_ds.R's own output directory) -> writes `leafcutter_results.tsv` (32,555 lines on the test set).
- `--majiq-dir` needs exactly one `*.deltapsi.tsv` in the directory -> writes `majiq_results.tsv` (126,949 lines). Not run against real MAJIQ output (licence-gated; MAJIQ itself is not installed here) — verified against Jutils' bundled `control_epileptic.deltapsi.tsv` instead, which is real MAJIQ output shipped with the tool.
- `--mntjulip-dir` needs `diff_spliced_groups.txt`, `diff_spliced_introns.txt` (DSR), `diff_introns.txt` (DSA), `intron_data.txt` -> writes `mntjulip_DSR_results_raw.tsv` / `mntjulip_DSA_results_raw.tsv`. The non-`_raw` `mntjulip_DS{R,A}_results.tsv` (with estimated per-sample PSI) is written **only** if the directory also has a `group_data.txt` with the extra estimated-PSI columns; the shipped test data has none, so only the `_raw` files appear. Pass whichever file exists to `heatmap`/`sashimi`/`venn-diagram`.
- `heatmap --tsv-file mntjulip_DSR_results_raw.tsv --meta-file data/mntjulip_meta_file.tsv --p-value 0.05 --q-value 1 --dpsi 0.2` and `sashimi ... --group-id g006855 --gtf data/gencode.vM17.annotation.clean.gtf.gz --shrink` both ran (PNGs written); `venn-diagram` combining all four tools' output TSVs also ran.

(Yang 2021 *Bioinformatics*) Useful when comparing multiple tools' outputs across publications or doing meta-analysis. The rMATS sashimi labels are per-sample junction counts and matched pysam.
