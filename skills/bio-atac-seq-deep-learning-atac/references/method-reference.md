# Sequence-Based Deep Learning for ATAC-seq: Method Reference

## Algorithmic Taxonomy

| Tool | Architecture | Training | Output | Strength | Fails when |
|------|-------------|----------|--------|----------|------------|
| chromBPNet (Pampari 2024 bioRxiv) | Bias-factorized CNN: a frozen Tn5 bias model plus an accessibility model; the bias model is trained on non-peak regions of the same library | Per cell type | Bias-corrected per-base profile + log total counts (heads: profile logits and log counts) | Purpose-built Tn5 bias factorization; pretrained ENCODE models exist | Too few reads/peaks; bias model that learned TF motifs (see QC gate). TF 2.8, CPU-only in the `chrombpnet` env |
| BPNet (Avsec 2021 Nat Genet 53:354) | Counts + profile dual-head CNN | TF ChIP-seq or ATAC | Per-base profile | Foundational; `bpnet-lite` reimplementation runs on PyTorch/GPU and also loads chromBPNet Keras `.h5` files | No built-in ATAC bias correction |
| scBasset (Yuan & Kelley 2022) | Basenji2-derived CNN with a learned per-cell embedding | Cells x peaks matrix | Per-cell, per-peak accessibility score; cell embedding | Sequence-informed scATAC embedding | Fixed architecture; TF <= 2.15 only (Keras 2); needs > ~20k peaks |
| Enformer (Avsec 2021 Nat Methods 18:1196) | Transformer, 196,608 bp input | Reference epigenomes (DNase/ATAC, ChIP, CAGE) | 896 bins x 128 bp x 5,313 human tracks | Distal regulation; pretrained | Fixed reference cell types; forward pass is GPU-heavy |
| Borzoi (Linder 2025 Nat Genet) | Enformer extension with RNA-seq | Multi-tissue | Sequence -> RNA + chromatin | Reported variant-effect benchmark on RNA | Not installed or executed here; no command shipped |
| DeepATAC / Basset (legacy) | Earlier CNNs | -- | Binary peak prediction | Historical context | Superseded; do not use for new work |
| tangermeme | PyTorch inference utilities | Any torch model | Ref/alt predictions, marginalization, DeepLIFT/SHAP | Batched GPU inference on torch models | Torch only: chromBPNet `.h5` must be loaded through bpnet-lite first |

Rules of thumb inherited from the original Skill (>= 50M deduplicated reads, >= 30k peaks, A100 ~24 h per cell type) are uncited and were not reproduced; the bounded CPU runs here (2.2k peaks, 2 epochs) prove workflow mechanics only, not model quality. Verify against the chromBPNet wiki, scBasset and Enformer/Borzoi papers before locking a pipeline.

## When Deep Learning Helps vs When Classical Pipelines Suffice

| Task | Classical | Deep Learning |
|------|-----------|---------------|
| Peak calling | MACS3 / Genrich (sufficient) | chromBPNet (overkill unless variant downstream) |
| Tn5 bias correction at TF motifs | TOBIAS ATACorrect (good) | chromBPNet (designed for low-complexity flanks and weak motifs) |
| Differential accessibility | DiffBind / DESeq2 (sufficient) | -- (no clear DL advantage) |
| GWAS variant effect prediction at causal SNPs | Limited (overlap heuristics) | chromBPNet / Enformer |
| Motif discovery from de novo data | MEME / HOMER (good) | chromBPNet + TF-MoDISco (can find composite/cooperative motifs) |
| Per-cell TF activity | chromVAR (sufficient at the cluster level) | scBasset (fine-grained cell states) |
| Cross-cell-type accessibility prediction | -- | Enformer (only option here) |

Classical pipelines remain primary for standard ATAC analysis. Deep learning enters when variant interpretation is the goal, prediction beyond observed data is needed, or bias-correction quality is paramount (low-input, FFPE, weak motifs).

## chromBPNet 1.0.1 pipeline contracts

[`scripts/chrombpnet_pipeline.sh`](../scripts/chrombpnet_pipeline.sh) holds the commands; these are the contracts it relies on.

- `prep splits`: `-tcr` = TEST chromosomes, `-vcr` = validation; all others train. `-tecr` does not exist.
- `prep nonpeaks -g -c -p -fl -o <prefix>` writes `<prefix>_negatives.bed` and refuses to run if `<prefix>_auxiliary/` exists. Nothing else creates the `-n` file that `bias train`/`train` need.
- Peaks, `pred_bw -r`, `contribs_bw -r` and `footprints -r` are 10-column narrowPeak/BED files (summit in column 10). A 3-column BED fails with `ValueError: cannot convert float NaN to integer`. MACS3 narrowPeak already qualifies; from a 3-column BED use `sort -k1,1 -k2,2n -u regions.bed | awk -v OFS='\t' '{print $1,$2,$3,".",0,".",0,0,0,int(($3-$2)/2)}'`. Duplicate windows make the bigWig writer fail with `entries ... out of order`, so keep one row per region, coordinate-sorted.
- `-b` is the bias threshold factor (0.5 ATAC, 0.8 DNase). Raise it and retrain the bias model if the bias model has negative peak correlation.
- `bias train` / `train` stop after training. `bias pipeline` / `pipeline` (and `qc`) append DeepSHAP interpretation, TF-MoDISco and PDF/HTML reports on up to 30,000 peaks; on CPU that tail dominates runtime (see Runtime). `--num-filters` does not exist; the flag is `-fil/--filters`; `-e/-es` cap epochs and early stopping.
- The bias model is trained on non-peak regions of the same library, not on a separate naked-DNA control. A bias model from another sample (ENCODE `model.bias_scaled.*.h5`) can substitute at degraded fidelity; the QC gate then matters more.
- Models: `model/models/chrombpnet_nobias.h5` (accessibility only, use for variants and motifs), `chrombpnet.h5` (with bias), `bias/models/bias.h5`.
- Outputs of `pred_bw -op <p>`: `<p>_chrombpnet_nobias.bw`, `<p>_bias.bw`, `<p>_chrombpnet.bw`, and with `-bw <observed.bw>` `<p>_<model>_metrics.json`. Only the requested regions carry signal.

### Runtime (fixture: 2,166 real GM12878 ATAC peaks re-cut into 5 pseudo-contigs, 2 epochs, 24-thread CPU shared with other jobs; not representative of production data)

| Step | Wall time |
|---|---|
| `prep splits` + `prep nonpeaks` | 8 s |
| `bias train` | 8.5 min |
| `train` (accessibility model) | 35 min |
| QC `pred_bw` (3 models, 403 test peaks) | 2.5 min |
| QC `footprints` (6 motifs, ~100 non-peak windows) | 12 min |
| Variant scoring, 30 SNPs | 3 min (forward + reverse-complement, 10 shuffles per SNP) |

The old `bias pipeline` and `pipeline` commands spent 61+ min in the DeepSHAP tail without exiting on this fixture (outputs were written first; the exit was never observed). Full-scale training is untested; expect training time to scale with peaks and epochs, and use `TRAIN_ARGS="-e 2 -es 1"` for a smoke run. A run cut short by a time limit can be continued with `RESUME=1`.

## QC gate

Run before any variant or motif use. Thresholds are the ones chromBPNet 1.0.1 prints in its own report (`make_html.py`, `marginal_footprinting.py`); [`scripts/chrombpnet_qc_gate.py`](../scripts/chrombpnet_qc_gate.py) applies them.

| Check | Source file | Pass |
|---|---|---|
| Bias model peak correlation | `model/evaluation/bias_metrics.json` `counts_metrics.peaks.pearsonr` | > -0.3 (train aborts below -0.5; -0.5..-0.3: inspect motifs for GC-rich hits) |
| Corrected model peak correlation | `pred_bw -bw` metrics on held-out peaks | > 0.5 |
| Tn5 bias response of the nobias model | `chrombpnet footprints` -> `*_max_bias_response.txt` | begins `corrected_` (max marginal footprint < 0.003 on all Tn5 motifs) |
| Bias model must not learn TF motifs | `chrombpnet bias qc` / `bias pipeline` TF-MoDISco of the bias model | no TF-like motifs (interpretation tail; run on suitable hardware) |

The 2-epoch fixture model correctly fails the correlation check (peak pearsonr 0.47 < 0.5); the gate reports that rather than passing it.

## In Silico Variant Effect Prediction

**Head semantics.** chromBPNet and BPNet return `(profile logits, log counts)`. The count fold change is `log2(exp(y_alt) / exp(y_ref)) = (y_alt - y_ref) / ln 2` computed on the LOG-count head. Taking `log2(y_alt / y_ref)` of the log-count output is wrong: on a planted-motif fixture (true log2FC -2.44) it returned -0.46 and the correct formula -2.26.

**Route A, variant-scorer (chromBPNet `.h5`, TF 2.8, CPU).**

```bash
git clone https://github.com/kundajelab/variant-scorer     # tested at 0e1e341; `chrombpnet snp_score` is not available
python variant-scorer/src/variant_scoring.py --model model/models/chrombpnet_nobias.h5 \
    --list variants.tsv --genome hg38.fa --chrom_sizes hg38.chrom.sizes --out_prefix variants --no_hdf5
```

- `variants.tsv`: 5 tab-separated, headerless columns `chr pos ref alt variant_id` (1-based pos), e.g. `chr1\t976669\tT\tC\tchr1:976669:T:C`. Four columns fail with `File has 4 columns but chrombpnet schema expects 5 columns`. `pybedtools` must be importable (install from conda, the pip build fails).
- Output `variants.variant_scores.tsv` (18 columns): `allele1_pred_counts`, `allele2_pred_counts`, `logfc` (= log2 alt/ref counts), `abs_logfc`, `jsd`, `logfc_x_jsd`, `abs_logfc_x_jsd`, and empirical p-values `*.pval` against a shuffled-sequence null.
- By default predictions are averaged over forward and reverse-complement input; `--forward_only` disables that. `-n` sets shuffled scores per SNP (default 10) and bounds p-value resolution.
- Calling variants: rank and threshold with the null-based `abs_logfc.pval` and `jsd.pval` (correct across variants, e.g. Benjamini-Hochberg), not with a fixed cutoff. `abs(logfc) > 1` is only a heuristic (2 of 200 real NA12878 SNPs in GM12878 peaks exceeded it on the ENCODE GM12878 model). Report effect size and significance together and state the model's QC status.

**Route B, PyTorch (GPU), any chromBPNet accessibility `.h5`.** `BPNet.from_chrombpnet(h5)` converts the Keras file without TensorFlow (outputs matched Keras to 4e-5 on 300 real windows). [`scripts/score_variants_torch.py`](../scripts/score_variants_torch.py) wraps it with `CountWrapper`, tangermeme `substitution_effect` and the correct formula; on 200 real SNPs it agrees with variant-scorer `--forward_only` to 1e-3 (max abs) and with the default forward+reverse-complement average at Pearson 0.96. Do not pass the raw `(profile, counts)` model to `substitution_effect`: it returns lists and `.sum` fails.

```python
from bpnetlite.bpnet import BPNet, CountWrapper
from tangermeme.variant_effect import substitution_effect
model = CountWrapper(BPNet.from_chrombpnet("chrombpnet_nobias.h5")).eval()   # log-count head only
y_ref, y_alt = substitution_effect(model, X, substitutions)   # X (N,4,2114); substitutions rows [example, position, new_base 0-3 ACGT]
log2fc = (y_alt - y_ref).squeeze() / np.log(2)
```

`tangermeme.marginalize.marginalize(model, X, motif)` inserts a motif (string, or one-hot `(1, 4, len)`) into background sequences and returns the predicted change. `marginal_predict` and `tangermeme.io.adapter` do not exist. `substitution_effect` returns a `PerturbationResult` that also unpacks as `y_ref, y_alt`. Marginal effect (ref vs alt at the SNP) answers the GWAS question; saturation ISM answers "which bases matter".

## Attributions and TF-MoDISco

The maintained package is `modisco-lite` (CLI `modisco motifs`, `modisco report`); the original `kundajelab/tfmodisco` `TfModiscoWorkflow` API is unmaintained and incompatible.

Two routes to modisco input (both executed):

- **Torch, GPU (recommended for speed):** [`scripts/attributions_to_modisco_npz.py`](../scripts/attributions_to_modisco_npz.py) writes `ohe.npz` and `attr.npz` (hypothetical DeepLIFT/SHAP contributions, layout `(N, 4, L)`; `(N, L, 4)` fails with `Window (500) cannot be longer than the sequences`). 1,500 real GM12878 peaks took 3-7 min on the GPU (busy shared host) (`deep_lift_shap` may warn when a convergence delta exceeds 0.001).
- **chromBPNet CLI, CPU:** `chrombpnet contribs_bw -m chrombpnet_nobias.h5 -r peaks.narrowPeak -g hg38.fa -c hg38.chrom.sizes -op contribs -pc counts profile` writes bigWigs plus `contribs.counts_scores.h5` / `contribs.profile_scores.h5` (`shap`, `projected_shap`, `raw`); feed them directly with `modisco motifs -i <h5>`. 12 peaks with both heads took 9 min on a loaded CPU (subsample peaks; the run only demonstrates the file handoff, 12 sequences give no motifs); there is no `shap_to_modisco` tool.

```bash
modisco motifs -s ohe.npz -a attr.npz -n 2000 -w 500 -o modisco_results.h5    # or: -i contribs.counts_scores.h5
modisco report -i modisco_results.h5 -o report/ -s report/ -m jaspar_or_hocomoco.meme -t   # needs `tomtom` (MEME suite) on PATH
modisco report -i modisco_results.h5 -o report/ -s report/ -m jaspar_or_hocomoco.meme -l   # tomtom-lite; no MEME needed, Euclidean distance
```

`-n` caps seqlets per metacluster; `-w` is the half-window around the peak centre (default 400); `-z` is the seqlet core size (default 20); `-s` on `report` is the image-link suffix; `-f` sets the seqlet flank. On the ENCODE GM12878 model with 1,500 peaks, `modisco motifs` (~5 min) and `report` (~1 min) found 14 patterns, 13 with a top JASPAR2024 hit at q < 0.05, including CTCF (q 2e-11), AP-1 (TGACTCA), ETS (CAGGAAGT), NF-kB/REL, IRF (GAAACTGAAAC), NFY (TGATTGG) and KLF/SP (GCCCCGCCC); validate discovered motifs against JASPAR/HOCOMOCO before publication. TF-MoDISco on 1M peaks was not tested.

## Enformer

Not in the Kipoi zoo; use `enformer-pytorch` (weights `EleutherAI/enformer-official-rough`, ~1.9 GB). Input is exactly 196,608 bp (token ids A,C,G,T,N = 0-4); output `(1, 896, 5313)`, 128 bp bins over the central 114,688 bp. Track names live in `targets_human.txt` (calico/basenji, `manuscripts/cross2020`); the model file has none (index 12 = `DNASE:GM12878`). [`scripts/enformer_variant_effect.py`](../scripts/enformer_variant_effect.py) scores SNPs on chosen tracks (forward strand only, pseudocount 1). Pad or trim to the fixed length; a cell type without a matching track needs an explicit proxy track that you document.

## scBasset

scBasset learns one embedding per cell (the projection kernel of shape `(32, n_cells)`) from the cells x peaks matrix; it does not need cluster aggregation, and pseudobulk data per cluster is a different workflow (chromBPNet). Environment `scbasset` (see SKILL.md), CPU. Validated on 10x PBMC 5k (26,041 peaks on chr1-4 x 4,609 cells): about 610 s per epoch on CPU; 3 epochs give only val AUC 0.53, so real use needs many more.

```bash
python scBasset/bin/scbasset_preprocess.py --ad_file data.h5ad --input_fasta genome.fa --out_path processed   # h5ad needs var columns chr,start,end
python scBasset/bin/scbasset_train.py --input_folder processed --out_path out --batch_size 128   # default 1000 epochs, early stopping (AUC, patience 50)
```

`scbasset_preprocess.py` crashes with `number sections must be larger than 0` unless validation and test each hold >= 1,000 peaks (roughly > 20k peaks in total; filter peaks by a low cell fraction such as 0.3% rather than 5% on small data). Keras 3 (TF 2.21 tested) fails at `ModelCheckpoint` (`filepath must end in .weights.h5`). The FASTA must contain every chromosome in the peak set.

## Reconciliation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| chromBPNet predicts a strong effect; MACS does not call a peak | Sequence model captures latent regulatory potential | Trust chromBPNet for variant effect; not for peak calling |
| Enformer differs from chromBPNet at the same locus | Different context (196 kb vs 2 kb), cell types, and read-outs (128 bp bins vs base resolution) | Report both with their context size |
| TF-MoDISco motifs differ from JASPAR | Sequence-derived vs ChIP-validated | Check JASPAR for confirmation; composites are expected |
| chromBPNet correction differs from TOBIAS ATACorrect | CNN vs k-mer bias model | Compare; TOBIAS remains standard for routine use |

**Operational rule:** for high-confidence variant prediction, agree across two approaches (chromBPNet plus Enformer); report a single-tool call as exploratory. On 200 real NA12878 SNPs the two read-outs (Enformer `DNASE:GM12878` vs ENCODE chromBPNet GM12878) had Spearman 0.50 (forward strand, 2 bins), which is why agreement must be checked per variant rather than assumed.

## References

- Pampari A et al 2024 bioRxiv 2024.12.25.630221 (chromBPNet; preprint)
- Avsec Z et al 2021 Nat Genet 53:354-366 (BPNet)
- Avsec Z et al 2021 Nat Methods 18:1196-1203 (Enformer)
- Linder J et al 2025 Nat Genet (Borzoi; consult the publication for volume/pages)
- Yuan H & Kelley DR 2022 Nat Methods 19:1088 (scBasset)
- Shrikumar A et al 2017 ICML (DeepLIFT)
- Schreiber J 2025 bioRxiv 2025.08.08.669296 (tangermeme)
- Shrikumar A et al 2018 arXiv:1811.00416 (TF-MoDISco)
- Kelley DR 2020 PLoS Comput Biol 16:e1008050 (Basenji2)
