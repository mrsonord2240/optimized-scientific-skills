# Tool Failure Modes, Reconciliation, and Common Errors

## Per-Tool Failure Modes

### rMATS-turbo: Junction Read Imbalance

**Trigger:** A cassette exon's flanking exons have unequal read mapping opportunity (very short upstream exon, repeat-overlapping flanks, or low-mappability regions).

**Mechanism:** rMATS' binomial model treats inclusion vs skipping junctions as having equal mappability. When mappability differs between the two junction types, the PSI estimate is biased.

**Symptom:** "Significant" rMATS calls with no concordant change in leafcutter or MAJIQ at the same locus; ΔPSI direction inconsistent with sashimi-plot intuition.

**Fix:** Run Shiba (Kubota 2025 *NAR*) which corrects junction-imbalance, or filter rMATS hits requiring concordant detection in leafcutter.

### SUPPA2: Sparse Empirical Null at Low Replicate Count

**Trigger:** n=2 vs n=2 (or n=3 vs n=2) design with `--method empirical`.

**Mechanism:** SUPPA2's empirical null is constructed from between-replicate ΔPSI distributions binned by transcript expression. With few replicates, the binned null is sparse and conservative-looking but actually under-calibrated.

**Symptom:** Inflated FDR (15-30% reported in benchmarks; not reproduced here); many "significant" hits don't replicate or validate.

**Fix:** Switch to leafcutter or Shiba for n<=3 designs; or use `--method classical` (Wilcoxon) for very low replicate count; reconcile against orthogonal tool.

### MAJIQ V3: Complex LSV Interpretation

**Trigger:** A gene with 4+ alternative splice sites at one node (e.g. one source, multiple acceptors).

**Mechanism:** A complete LSV at a single node lists all observed junctions; PSI is per-junction within the LSV, not "PSI of one event."

**Symptom:** Reporting "PSI of the gene" doesn't make sense; per-junction PSIs sum to 1 across the LSV but no single number represents the gene.

**Fix:** Use VOILA to visualize the LSV graph and identify which junction(s) shifted; for cassette-style reporting, derive equivalent PSI from sum of inclusion-junctions / total junctions in the LSV.

### leafcutter: Cluster Topology Arbitrariness

**Trigger:** A cluster has 4+ introns sharing splice sites with non-canonical topology (e.g. mixed cassette + alternative donor + IR).

**Mechanism:** leafcutter clusters introns by shared splice sites; complex topologies don't map onto SE/A5SS/A3SS taxonomy and cluster-level "ΔPSI" hides which intron drove the change.

**Symptom:** Significant cluster-level p-value but multiple introns showing different effect-size directions.

**Fix:** Inspect the cluster in leafviz; report per-intron effect sizes (`effect_sizes.txt`); for canonical event reporting, map to SE/A5SS/A3SS via flanking exon coordinates manually.

## Reconciliation: When rMATS and leafcutter Disagree

The two most common short-read tools answer slightly different questions: rMATS classifies on annotated event templates; leafcutter classifies on observed cluster usage. Disagreement is informative.

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| rMATS sig, leafcutter not sig | rMATS junction read imbalance OR rMATS event hits an annotation that leafcutter clustered differently | Inspect locus in IGV; check Shiba |
| leafcutter sig, rMATS not sig | Novel junction not in rMATS annotation; rMATS `--novelSS` may have missed it | Check `--novelSS` was on; rerun if not |
| Both sig, opposite ΔPSI direction | Event class mismatch (e.g. rMATS calls SE positive but leafcutter sees A5SS shift in same cluster) | Manually map cluster topology to event class |
| Both sig, same direction | High-confidence call | Report; cross-validate with sashimi-plot |

**Operational rule:** for high-confidence reporting, require concordant detection in two tools from different algorithmic families (event-based + cluster-based, or LSV + isoform-based).

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `error: GTF gene_id parsing` (rMATS) | rMATS expects GENCODE-style gene_id; some Ensembl GTFs use different attribute order | `gffread input.gff3 -T -o standardized.gtf` |
| `KeyError: 'IJC_SAMPLE_1'` (rMATS parsing) | Output column missing; sometimes occurs when --statoff combined with novel events on older versions | Update rMATS-turbo to >=4.3.x; re-run |
| `MAJIQ: too few reads at junction` | Default `--minreads 10 --minpos 3` filters out the locus | Lower thresholds for low-coverage data; document filtering |
| `leafcutter: dispersion estimation failed` | Cluster has all-zero counts in one group | Pre-filter with leafcutter_ds.R `--min_samples_per_group 3 --min_samples_per_intron 5` |
| `SUPPA2: empirical p computed on N=4 nulls` | Insufficient replicates for empirical mode | Switch `--method classical` (Wilcoxon) for very low replicate count |
| `regtools: invalid CIGAR` | Non-BAM-spec read in input | Filter with `samtools view -h -F 0x100 -F 0x800` (drop secondary/supplementary) |
| `STAR: too many SJs` (in pass 2) | Cohort SJ.out.tab too large | Filter to junctions seen in >=3 samples or with >=3 unique reads before merging |

