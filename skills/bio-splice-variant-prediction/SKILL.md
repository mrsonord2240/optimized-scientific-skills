---
name: bio-splice-variant-prediction
description: Predicts whether a DNA variant alters mRNA splicing using sequence-based deep-learning tools — SpliceAI (10kb context dilated CNN, clinical default), Pangolin (multi-tissue), MMSplice (modular per-region CNN with calibrated ΔPSI), SpliceTransformer (tissue-aware transformer), CI-SpliceAI (SpliceAI retrained on all isoforms), SpliceVault (empirical 300K-RNA lookup of likely mis-splicing outcomes), CADD-Splice (composite score). Applies the ClinGen SVI 2023 framework as research-use computational evidence for ACMG/AMP interpretation (PVS1, PP3, BP4 codes), HGVS splicing nomenclature (c.123+1G>A, c.123-3T>G, r.spl?), extended-window scoring for deep-intronic pseudoexons, tissue-specific predictions, branchpoint-variant caveats, and a splice-switching ASO design checklist. Use when interpreting splice impact of variants, prioritizing VUS, identifying deep-intronic candidates, or planning ASOs.
tool_type: python
primary_tool: SpliceAI
license: MIT
---

## Version Compatibility

Checked 2026-09-20 on: SpliceAI 1.3.1 (TensorFlow 2.21, setuptools 80.x), Pangolin 1.0.2 (torch 2.13, PyVCF3 1.0.0, gffutils 0.14, pyfastx 2.3.1), MMSplice 2.4.0, CI-SpliceAI 1.2.2 (TensorFlow 2.15, keras 2.15), SpliceTransformer (GitHub `main`, torch 2.11, sinkhorn-transformer 0.11.4 + axial-positional-embedding 0.2.1), pysam 0.24.1, pandas 2.3, GENCODE v45.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Splice Variant Prediction

Predict whether a DNA variant alters mRNA splicing. **Distinct from "variant pathogenicity" generally**: a variant can be a strong splice disruptor without being pathogenic for the gene's standard mechanism, or pathogenic for reasons orthogonal to splicing. Splice prediction asks specifically: does this variant change splice-site usage?

**Scope.** Everything here is research-use decision support: computational evidence for an expert to weigh, not a diagnosis, a treatment recommendation, or a substitute for RNA validation and clinical-laboratory sign-off. Never write a prediction up as a patient result.

## Predictor Taxonomy

| Family | Architecture | Output | Fails when |
|--------|--------------|--------|------------|
| Context-aware CNN | 10 kb dilated ResNet | Per-position donor/acceptor probability | Long-range (>5 kb) regulatory effects; tissue-specific events |
| Tissue-aware CNN/transformer | Same arch + multi-tissue training | Per-tissue ΔPSI | Tissue not in training set; novel cell types |
| Modular per-region CNN | Separate sub-models for 5'ss/3'ss/exon/intron | Calibrated quantitative ΔPSI | Atypical events; complex multi-junction effects |
| Foundation transformer | Pretrained on broad genomic context | Splice probability or ΔPSI | New tools; less battle-tested |
| Empirical lookup | Public RNA-seq event database | Top-N most likely mis-splicing outcomes | Variant types not represented in training cohorts |
| Composite score | Blend of multiple predictors | Single scaled score | When component predictors disagree internally |

## Tool Selection Matrix

| Tool | Best for | Output | When to use | Fails when |
|------|----------|--------|-------------|------------|
| SpliceAI | Clinical screening; canonical splice site disruption | Delta score 0-1 | Default for ACMG variant classification | Tissue-specific events; effect site farther than -D from the variant |
| Pangolin | Tissue-aware predictions | Max gain/loss over 4 tissue models (per-tissue via Python) | When disease tissue is known (brain, heart, liver, testis) | Tissue not in 4-tissue training set; `-m True` (see Pangolin) |
| MMSplice | Quantitative ΔPSI | Δlogit_psi | Research where calibrated effect-size matters | Atypical events outside cassette-exon model; variants far from an exon |
| SpliceTransformer | 2024+ benchmark improvements | Splice score + 15 tissue-specificity flags (SNVs) | When transformer foundation models outperform CNN on benchmark variant sets | New (2024); limited clinical adoption |
| CI-SpliceAI | Second SpliceAI-family opinion | SpliceAI-style delta scores | Extended windows; cross-check | Did not flag GLA c.639+919G>A here (0.05) |
| SpliceVault | Empirical mis-splicing outcome | Top-N events at the affected splice site | Predicting consequence (skip vs cryptic) of canonical-disrupting variants | Variants not at catalogued splice sites |
| CADD-Splice | Single composite score | Scaled C-score (PHRED) | Clinical pipelines wanting one number | When knowing which sub-component drove the score is needed |

Methodology evolves; verify benchmarks (Smith & Kitzman 2023 *Genome Biol* 24:294; You et al 2024 *Nat Commun*) and ClinGen SVI splicing recommendations before reporting any interpretation. Concordance across SpliceAI + Pangolin + MMSplice is the strongest computational evidence; discordance flags need RNA validation.

## Decision Tree by Use Case

| Use case | Recommended approach |
|----------|----------------------|
| Single variant, ACMG-style evidence | SpliceAI default 50nt + ClinGen SVI 2023 thresholds |
| Tissue-specific question (brain disease, cardiomyopathy) | SpliceAI + Pangolin per-tissue (`pangolin_tissue.py`) or SpliceTransformer |
| Unsolved Mendelian case (suspect deep-intronic) | SpliceAI `-D 500` on every intronic candidate, both pseudoexon ends reported; SpliceVault for the consequence |
| VUS panel screening | SpliceAI + Pangolin + MMSplice concordance (`splice_parsers.py`) |
| Predict consequence of canonical-disrupting variant | SpliceVault top-N empirical events |
| Branchpoint variant suspected | Predictors are weak here; RNA validation (see Branchpoint Variant Detection) |
| Splice-switching ASO design | Checklist only (see ASO section); no executable here |
| Validate predicted splice change in patient | RNA-seq + FRASER2 (see outlier-splicing-detection) |
| Pseudoexon prediction in deep intron | SpliceAI extended window; CI-SpliceAI as second opinion; require RNA validation |

## ClinGen SVI 2023 Framework

The ClinGen Sequence Variant Interpretation (SVI) splicing subgroup (Walker 2023 *Am J Hum Genet*) extended the ACMG/AMP 2015 framework with explicit splice-prediction rules.

| Evidence code | Threshold | Notes |
|----------------|-----------|-------|
| **PP3** (supporting pathogenic) | SpliceAI delta >= 0.20 | ClinGen SVI: apply at **supporting** weight (not standalone) |
| **BP4** (supporting benign) | SpliceAI delta <= 0.10 | ClinGen SVI: apply at **supporting** weight |
| **PVS1** (very strong null) | Canonical +/-1, +/-2 site disruption with predicted LoF + NMD | Requires gene where LoF is established mechanism (Abou Tayoun 2018 *Hum Mutat* PVS1 decision tree) |
| **PS3 / BS3** (functional) | RNA evidence (RT-PCR, RNA-seq, minigene) | Supersedes computational evidence |
| *(no code)* | delta 0.10-0.20, or **no score** (`.`, no `SpliceAI=` tag) | Inconclusive / `not_scored`: a missing score is never BP4 |

Boundaries are inclusive as published: 0.10 -> BP4, 0.20 -> PP3. SpliceAI scores have two decimals, so exact 0.20/0.50/0.80 values are common.

**Operational rules:** Computational evidence (PP3/BP4) is *supporting*, not standalone. Higher SpliceAI cutoffs (0.5, 0.8) increase precision but are the tool's own tiers (Jaganathan 2019), NOT ClinGen-endorsed evidence-strength upgrades — reaching moderate/strong requires functional/RNA evidence (PS3/BS3), not a higher SpliceAI score alone. Splicing variants benefit from concordance across SpliceAI + Pangolin + MMSplice. RNA validation supersedes prediction. Always log SpliceAI version, distance window, and reference transcript. SpliceAI alone is **not sufficient** for PVS1 (canonical site disruption also needs gene-level LoF context) or for non-canonical positions.

## Install (Linux/WSL; checked 2026-09-20)

```bash
pip install spliceai tensorflow "setuptools<81"     # SpliceAI imports pkg_resources; weights and grch37/grch38 gene tables ship in the package

# Pangolin: GitHub, NOT the PyPI/bioconda 'pangolin' (SARS-CoV-2 lineages). Its setup.py declares no dependencies:
pip install torch gffutils pyfaidx pyfastx biopython pandas "PyVCF3==1.0.0"   # PyVCF3 >= 1.0.2 breaks vcf.parser._Info
git clone https://github.com/tkzeng/Pangolin && pip install ./Pangolin     # also puts create_db.py on PATH

pip install mmsplice       # needs setuptools<81 and a cyvcf2 built for the installed numpy (else: pip install --force-reinstall --no-deps cyvcf2)
```

Reference: an **upper-case** FASTA of the same build as the VCF (GENCODE `GRCh38.primary_assembly.genome.fa` is upper-case; UCSC/Ensembl soft-masked FASTAs are not, see Common Errors) and a GENCODE **GTF** (not GFF3) for Pangolin/MMSplice. SpliceVault, SpliceTransformer, CI-SpliceAI and CADD installs are in their sections.

**Data handling.** SpliceVault (remote tabix), CADD (API) and VariantValidator/Mutalyzer send variant coordinates to third-party servers; use local files or skip them for identifiable patient data.

## SpliceAI Workflow

**Goal:** Annotate VCF variants with per-variant delta scores for splice-site change.

**Approach:** Run `spliceai` CLI with reference genome and annotation; parse the INFO field. **SpliceAI is human-only** (`-A grch37` or `-A grch38`, which must match the FASTA); the model was trained on GENCODE human and does not directly transfer to mouse, fly, or other species. For mouse, retrained variants exist (e.g. mouseSpliceAI); for other species, use Pangolin (4 species: human, mouse, rat, rhesus macaque) or accept that prediction will be unreliable.

```bash
spliceai \
    -I input.vcf \
    -O output.vcf \
    -R GRCh38.primary_assembly.genome.fa \
    -A grch38 \
    -D 50 \
    -M 0
```

`-D` = **maximum distance between the variant and the gained/lost splice site** (default 50). `-M 0` (default) returns raw scores; `-M 1` masks splice gains at annotated sites and losses at unannotated sites (annotation = SpliceAI's canonical GENCODE table). Output INFO format: `SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`. Delta score = max(DS_AG, DS_AL, DS_DG, DS_DL). DS labels: AG = acceptor gain, AL = acceptor loss, DG = donor gain, DL = donor loss; DP = position of that site relative to the variant. One record can carry several annotations (readthrough/overlapping genes such as RPL36A-HNRNPH2 add rows): reduce per variant, and choose the MANE gene when genes disagree.

**SpliceAI exits 0 on records it cannot score.** REF mismatch, a deletion longer than 2*D, or no overlapping gene give a record with **no `SpliceAI=` tag** (stderr warning only); N-containing sequence gives `SpliceAI=...|.|.|.|.|.|.|.|.` (all dots). A `<DEL>` symbolic ALT crashes it. Assert REF against the FASTA first and compare output to input: a failed variant must not read as benign or vanish.

Parse, label and check with `examples/splice_parsers.py` (`.` -> NaN -> `not_scored`, boundaries inclusive, unmatched records reported):

```python
from splice_parsers import parse_spliceai_vcf, classify_delta, read_input_vcf, unscored_report

s = parse_spliceai_vcf('output.vcf')                       # one row per variant x ALT x gene
s['acmg_evidence'] = classify_delta(s['delta_max'])        # BP4 / inconclusive / PP3_supporting[_prec0.5|_prec0.8] / not_scored
print(unscored_report('input.vcf', s['key'], 'SpliceAI'))  # input records with no SpliceAI= tag
```

`examples/spliceai_clingen_classify.py input.vcf genome.fa --build grch38` runs the whole flow (default window, wider re-run of everything scoring < 0.20 including 0.00, PP3/BP4 labels, skipped-record report). `python test_splice_parsers.py` in `examples/` self-tests the parsers and boundary labels on real SpliceAI/Pangolin/MMSplice output (`test_data/`, GRCh38 panel: TP53 c.673-2A>G, DMD c.9563+1G>A / c.31+1G>A, OTC c.386+5G>A, GLA c.370-1G>A, GLA c.639+919G>A, three ClinVar-benign GLA variants). Expected SpliceAI D50 delta_max on that panel: TP53 1.00, DMD c.9563+1G>A 0.99, DMD c.31+1G>A 0.95, OTC c.386+5G>A 0.94, GLA c.370-1G>A 1.00, GLA c.639+919G>A 0.30, the three benign GLA variants 0.00.

## Pangolin for Tissue-Specific Prediction

**Goal:** Get tissue-aware splice impact predictions when disease tissue is known.

**Approach:** Run Pangolin CLI with VCF + reference + gffutils annotation database built from a GENCODE **GTF**.

```bash
# annotation DB: Pangolin's own script; default keeps only Ensembl_canonical transcripts (several minutes for the full GTF)
create_db.py gencode.v45.annotation.gtf          # -> gencode.v45.annotation.db   (add --filter None to keep every transcript)

pangolin input.vcf GRCh38.primary_assembly.genome.fa gencode.v45.annotation.db pangolin_output -d 50 -m False
pangolin input.vcf GRCh38.primary_assembly.genome.fa gencode.v45.annotation.db pangolin_d500 -d 500 -m False -s 0.2   # every site with |change| >= 0.2
```

Do not build the DB with `gffutils.create_db('gencode.v45.annotation.gff3', ...)`: on real GENCODE GFF3 it raises `ValueError: Duplicate ID CDS:ENST...`. Only substitutions and simple indels are scored (REF or ALT one base); a multi-allelic record scores **ALT[0] only, with no warning** (`bcftools norm -m-` first); output uses Ensembl gene IDs (`gene|pos:gain|pos:loss|Warnings:`, versioned IDs from GENCODE) and Pangolin compares REF case-sensitively.

**Tissue-specific output.** The CLI prints the maximum gain and maximum loss over the four tissue models (heart, liver, brain, testis), not per-tissue scores. Per-tissue values come from `examples/pangolin_tissue.py`, which reuses the installed models:

```bash
python pangolin_tissue.py chr17 7674292 T C --strand - --fasta GRCh38.primary_assembly.genome.fa
# heart: gain +0.38 at +47, loss -0.86 at -2 | liver ... | brain: gain +0.38, loss -0.87 | testis: gain +0.72, loss -0.90
# tissue maximum (what the CLI reports): gain +0.72, loss -0.90
```

The models were trained on **4 tissues: brain, heart, liver, testis** (Zeng & Li 2022 *Genome Biol*) and extrapolate poorly to tissues outside this set: for those, fall back to SpliceAI.

**`-m True` is not safe for variant screening** (README default is True; measured on GENCODE v45, Pangolin 1.0.2, `-d 50`; losses are Pangolin scores, negative):

| Variant | `-m False` | `-m True`, canonical-only DB | `-m True`, all-transcript DB |
|---------|-----------|------------------------------|------------------------------|
| GLA c.639+919G>A (pseudoexon gain, GRCh38 X:101399747 C>T) | +0.23 at -3 (and +0.23 at +53 with `-s 0.2`) | +0.23 | **0.00** (the pseudoexon is annotated in NMD transcript ENST00000493905) |
| 5 DMD splice sites annotated only in non-canonical transcripts (losses) | -0.38 to -0.81 | **0.00 (all five erased)** | -0.38 to -0.81 |
| OTC c.386+5G>A (canonical donor loss; OTC overlaps ENSG00000250349) | -0.72 at -5 | **0.00 (erased)** | **0.00 (erased)** |
| TP53 c.673-2A>G (canonical acceptor loss) | -0.90 at -2 | -0.90 | -0.90 |

The masking zeroes gains at every annotated site and losses at every unannotated one, so neither DB is neutral, and with overlapping genes an earlier gene's mask is applied in place to the arrays the next gene uses (an OTC-only DB keeps -0.72). Run `-m False` and read gains/losses against the annotation yourself; if you also run `-m True`, any loss present at `False` and absent at `True` needs review, not trust. Pangolin `-m False` with the canonical DB agreed with SpliceAI on all nine panel variants (loss/gain >= 0.2 at the same sites).

## SpliceVault for Empirical Mis-Splicing Outcomes

**Goal:** Predict the *type* of mis-splicing (exon skipping vs cryptic site activation) given a canonical-disrupting variant.

**Approach:** Query the public SpliceVault table (GRCh38 SNVs at annotated splice sites, one row per transcript) by remote tabix; `examples/splicevault_lookup.py` returns the parsed Top-4 events. Needs `pysam` and `pandas`; downloading the 884 MB `SpliceVault_data_GRCh38.tsv.gz` is optional (`--tsv local.gz`). Bulk SQL dumps at `storage.googleapis.com/misspl-db-data` are requester-pays (billing project required); the web portal is `kidsneuro.shinyapps.io/splicevault`; the same table backs the SpliceVault VEP plugin (Ensembl 111+).

```bash
curl -O https://ftp.ensembl.org/pub/current_variation/SpliceVault/SpliceVault_data_GRCh38.tsv.gz.tbi   # index only; rows stream over HTTP
python splicevault_lookup.py chr17 7674292 T C --transcript ENST00000269305 --tbi SpliceVault_data_GRCh38.tsv.gz.tbi
# ENST00000269305 Acceptor_loss at chr17:7674291; SpliceAI delta 1; out-of-frame Frameshift:3/4; 199336 samples
# Top1 CA +47 0.4% Frameshift | Top2 CA -50 0.08% Frameshift | Top3 ES 7 0.03% Frameshift | Top4 CA -70 0.03% inFrame
```

Event types: ES = exon skipping (impact = skipped exon number), CA/CD = cryptic acceptor/donor (impact = nt from the annotated site, transcript sense). Checked against SpliceAI on the same variant: Top1 `CA +47` is SpliceAI's acceptor gain (DS_AG 0.85, DP_AG +47). `out_of_frame` is the fraction of the Top-N events that shift the frame; >= 3/4 in-frame suggests a variant that may be splice-rescued rather than LoF. An empty result means the variant is not a catalogued splice-site SNV (or REF/ALT/build is wrong).

SpliceVault (Dawes 2023 *Nat Genet*; 335,663 RNA-seq samples) reports that the **300K-RNA Top-4 events** predict variant-associated mis-splicing with 92% sensitivity, identifying 96% of exon-skipping events and 86% of cryptic-site events in 140 clinical cases with RNA testing. Use it when the question is not "will splicing change?" but "what specific aberrant splicing will occur?".

## MMSplice for Calibrated ΔPSI

**Goal:** Predict quantitative ΔPSI (not just probability of disruption) for cassette exons.

**Approach:** Score variant impact on each splicing region (5'ss, 3'ss, exon, intron-3'/5') and combine.

```python
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save

dl = SplicingVCFDataloader(
    gtf='gencode.v45.basic.annotation.gtf',
    fasta_file='GRCh38.primary_assembly.genome.fa',
    vcf_file='input.vcf'
)

model = MMSplice()
predict_save(model, dl, 'mmsplice_predictions.csv', pathogenicity=True)
```

MMSplice (Cheng 2019 *Genome Biol*) reports Δlogit_psi per variant. Useful when calibrated effect sizes matter (research) more than probability of disruption (clinical screening). The CSV has **one row per variant x exon x transcript**, ID `chrom:pos:ref>alt`, and no row for a variant that is not near an annotated exon (the deep-intronic GLA c.639+919G>A got none, silently): reduce with `parse_mmsplice_csv()` (largest |delta_logit_psi| per variant) and report variants missing from the output. A truncated or mismatched FASTA yields empty sequences and assertion errors (`Input sequence ... intron length cannot be longer than the input sequence`). Companion **MTSplice** (Cheng 2021 *Genome Biol*) adds tissue-specific Δψ predictions.

## SpliceTransformer, CI-SpliceAI and CADD-Splice

**SpliceTransformer** (You 2024, GitHub `ShenLab-Genomics/SpliceTransformer`; SNVs only, indels skipped silently; hg38 or hg19). Weights are not in the repo: a public Google Drive file (`gdown 1d8n4vHDSbXqpPc_JFEswLomSUDBgHvno`, 120 MB) saved as `model/weights/SpTransformer_pytorch.ckpt`.

```bash
pip install torch sinkhorn-transformer "axial-positional-embedding==0.2.1" "PyVCF3==1.0.0" pyensembl gffutils pyfaidx pandas tqdm gdown
# axial-positional-embedding 0.2.1 is required: newer versions rename pos_emb weights and load_state_dict fails
# put a chr-prefixed hg38.fa and GENCODE hg38.annotation.gtf.gz in data/data_package/, then index the GTF once
# (the annotation name is hard-coded in sptransformer.py, whatever GENCODE release the GTF is):
pyensembl install --reference-name hg38 --annotation-name gencode.v38 --gtf data/data_package/hg38.annotation.gtf.gz
python sptransformer.py -I input.vcf -O out.csv --reference hg38
```

Output: `score` (largest change in splice-site probability, 0-1) plus Y/N flags for 15 tissues (tissue-specific usage change above the paper's 95% threshold). Panel result: TP53 c.673-2A>G 1.00, DMD c.9563+1G>A 0.99, GLA c.370-1G>A 0.98, GLA c.639+919G>A 0.38, ClinVar-benign GLA 0.01-0.07. Runtime is ~25 s per variant on CPU (3 min 52 s for 9 variants) and 12 s for the same 9 on a GPU: the RTX 5070 Ti (Blackwell, sm_120) needs a CUDA 12.8 torch build (`pip install torch --index-url https://download.pytorch.org/whl/cu128`, checked with torch 2.11.0+cu128); scores were identical on CPU and GPU.

**CI-SpliceAI** (Strauch 2022; SpliceAI retrained on all GENCODE isoforms). `pip install "tensorflow-cpu==2.15.*" "keras<3" cispliceai` (keras 3 is not supported; 1.2.2 checked), then:

```bash
cis-vcf -a grch38 -d 500 --all -i input.vcf -o ci_output.vcf GRCh38.primary_assembly.genome.fa
```

Output INFO tag is `CISpliceAI=ALLELE|ENSEMBL_GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_...` (SpliceAI layout; parse it by replacing `SpliceAI=`); `-a grch37` for GRCh37. On the GRCh37 panel it flagged every canonical variant (delta_max 0.64-1.00; SpliceAI 0.95-1.00) but scored the GLA c.639+919G>A pseudoexon 0.05 at `-d 500` where SpliceAI gave 0.30, so it is a cross-check, not a rescue.

**CADD-Splice** is the standard CADD score from v1.6 on (splice predictors are inputs to the model). The public API needs no account:

```bash
curl -s "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/17:7674292_T_C"
# [{"Alt":"C","Chrom":"17","PHRED":"34","Pos":"7674292","RawScore":"6.512109","Ref":"T"}]
```

PHRED on the panel (GRCh37-v1.7): canonical splice variants 33-35 (DMD c.31+1G>A 33, DMD c.9563+1G>A 34, GLA c.370-1G>A 35), GLA c.639+919G>A 14.9, benign GLA rs2071228 10.7. It does not say which component drove the score; report it beside, not instead of, SpliceAI/MMSplice.

## HGVS Splicing Nomenclature

Following den Dunnen 2016 *Hum Mutat*:

| Notation | Meaning |
|----------|---------|
| `c.123+1G>A` | +1 of intron downstream of exon ending at cDNA position 123 (canonical 5'ss G) |
| `c.123+5G>A` | +5 position of donor (consensus region) |
| `c.124-1G>A` | -1 of acceptor (canonical AG) |
| `c.124-3T>G` | -3 of acceptor (Py-tract / BPS region) |
| `c.124-50A>G` | Deep-intronic; may activate cryptic site |
| `r.123_456del` | RNA-level deletion (predicted exon skipping) |
| `r.spl?` | Unknown splice consequence |
| `r.0?` | No detectable RNA |
| `p.0?` | Unknown protein consequence |
| `p.(=)` | No predicted protein change (silent) |

Validation tools: VariantValidator (Freeman 2018 *Hum Mutat*), Mutalyzer 2 (Lefter et al 2021 *Bioinformatics* 37:2811-2817). Check REF against the reference FASTA before scoring; a wrong REF is skipped without failing the run.

## Extended-Window Scoring for Deep-Intronic Variants

`-D` is the distance from the variant to the gained/lost site, not the distance to the nearest canonical site. A pseudoexon-creating variant usually scores **at the variant itself** (it creates a donor or acceptor there); a wider `-D` adds the *other* pseudoexon boundary. Measured on GLA c.639+919G>A (919 nt from exon 4): `-D 50` DS_DG 0.30 at the variant; `-D 500` and `-D 2000` add DS_AG 0.22 at +53 (the far end), delta_max unchanged at 0.30. Pangolin `-d 500 -s 0.2 -m False` lists both ends (+0.23 at -3 and +0.23 at +53).

| Window | Tradeoff |
|--------|----------|
| -D 50 (default) | Fast; captures canonical-site disruption and gains created at the variant |
| -D 500 | Adds far-end boundaries and distant sites (e.g. DMD c.31+1G>A donor gain 0.08 at -10 with -D 50 becomes 0.40 at +69) |
| -D 2000 | Maximum sensitivity; more distant, weaker sites |

For unsolved cases run `-D 500` on every intronic candidate that scores below 0.20 (including 0.00), and report DS **and DP of both ends** of any candidate pseudoexon:

```bash
spliceai -I candidates.vcf -O output_D500.vcf -R genome.fa -A grch38 -D 500 -M 0
```

Pseudoexon creation in deep introns explains a substantial fraction of unsolved Mendelian disease alleles in current cohorts (estimates 5-15% across studies; specific quantitative range will vary by cohort and panel — verify against current literature). Disease examples: CFTR 3849+10kbC>T, USH2A c.7595-2144A>G, CEP290 c.2991+1655A>G (LCA10), GLA c.639+919G>A (Fabry). CI-SpliceAI (see above) is a second opinion, not a rescue.

## Concordance Across Predictors

`build_concordance()` in `examples/splice_parsers.py` parses the three tools' real outputs, joins them on a normalised `chrom:pos:ref>alt` key (`chr` prefix dropped), keeps every input variant (outer join: a variant a tool skipped shows NaN and lowers `n_scored`) and never drops a tool's missing output silently:

```python
from splice_parsers import build_concordance
T = build_concordance('input.vcf', spliceai_vcf='output.vcf', pangolin_vcf='pangolin_output.vcf',
                      mmsplice_csv='mmsplice_predictions.csv')
print(T[['spliceai_delta', 'pangolin_score', 'delta_logit_psi', 'n_scored', 'n_above', 'concordance']])
```

Thresholds per tool: SpliceAI delta >= 0.2, Pangolin |score| >= 0.2, MMSplice |delta_logit_psi| >= 1.0 (Skill conventions for Pangolin/MMSplice; only the SpliceAI one is ClinGen). Multi-allelic records must be split first (`bcftools norm -m-`): Pangolin scores ALT[0] only.

| Concordance label | Meaning | Action |
|-------------------|---------|--------|
| `all_predict_disruption` | every tool that scored is above threshold | PP3 (supporting); strong candidate for RNA validation (PS3) |
| `majority_predict_disruption` | more than half above (e.g. 2/3) | PP3 (supporting) |
| `discordant` | some, but not most, above | Report inconclusive; flag for RNA validation |
| `none_predict_disruption` | all below | BP4 (supporting) only if SpliceAI <= 0.10 |
| `insufficient_tools` | fewer than 2 tools scored | Find out why (skipped records) before interpreting |

Discordance is the most informative pattern — variants where one model sees impact and others don't are high priority for RNA validation. Check `n_scored`: on the test panel the GLA c.639+919G>A pseudoexon is `all_predict_disruption` from only two tools because MMSplice returned nothing.

## Branchpoint Variant Detection

All current tools are **weak at branchpoint variants** because the BPS motif (yUnAy) has low information content; SpliceAI captures only some. Published branchpoint-specific methods: BPP (Zhang 2017 *Bioinformatics* 33:3166), LaBranchoR (Paggi & Bejerano 2018 *RNA* 24:1647), SVM-BPfinder (Corvelo 2010 *PLoS Comput Biol*), BPHunter (Zhang 2022 *PNAS*; web server + standalone). **None was run here**: BPHunter's reference datasets were announced at `hgidsoft.rockefeller.edu/BPHunter/standalone.html`, which redirects to a GitHub page that returns 404 (2026-09-20), so its standalone script cannot be run. Recommendation: when SpliceAI delta is borderline (0.1-0.3) for a variant in the BPS region (-18 to -40 from 3'ss), treat all predictors as uninformative and require RNA validation; BPHunter's web server can be tried as a supplement.

## Splice-Switching ASO Design

**Goal:** Plan antisense oligonucleotides to modulate splicing therapeutically (e.g. SMA ISS-N1, DMD exon skipping). This is a **checklist for discussion with a design platform, not an executable pipeline**; no design tool is provided or was run here, and SpliceAI cannot model oligonucleotide occlusion (on N-masked sequence it returns `.` for every score).

1. Identify the target ESE/ESS/ISE/ISS region (splice-site strength with MaxEntScan; motif databases)
2. Design candidate 18-22 nt ASOs spanning the regulatory element
3. Filter for RNA accessibility (avoid stable hairpins) with RNAfold
4. Transcriptome-wide sequence search for off-target binding (<=16/20 nt match to any non-target transcript)
5. Avoid TLR9 immunostimulatory CpG motifs
6. Score the *variant/sequence change* you aim to correct or mimic with SpliceAI, never the oligo

Chemistry choices: 2'-MOE-PS (nusinersen-like, CNS, intrathecal); PMO (DMD ASOs, systemic IV); GalNAc-conjugated (hepatic targeting). Approved precedents: **nusinersen** (SMA ISS-N1 occlusion, exon 7 inclusion); **risdiplam** (small-molecule SMN2 splicing modulator); **eteplirsen/golodirsen/casimersen/viltolarsen** (DMD exon skipping). Design references: Hua 2008 *AJHG*; Roberts et al 2023 *Nat Rev Drug Discov* 22:917 (DMD therapeutic approaches).

## Per-Tool Failure Modes

### SpliceAI: Effect Site Farther Than -D

**Trigger:** The gained/lost splice site lies more than `-D` nt from the variant (deep-intronic pseudoexon whose creating change is far from its boundaries, or a distant cryptic site).

**Mechanism:** `-D` limits the variant-to-site distance; nothing outside is scored.

**Symptom:** Known pathogenic deep-intronic variant scores low (<0.2) or 0.00; no pseudoexon detected.

**Fix:** Re-run with `-D 500` (or `2000`) and read the DP of every site; a variant that already scores at itself (like GLA c.639+919G>A, 0.30) is found at `-D 50`. CI-SpliceAI is a second opinion only.

### SpliceAI: Tissue Agnosticism

**Trigger:** Variant in a tissue-specific gene (NEFM in neurons, MAPT brain, DMD muscle isoforms).

**Mechanism:** SpliceAI is trained on aggregate GENCODE annotation; tissue-specific events with weak constitutive use score low.

**Symptom:** Tissue-specific pathogenic variant has low SpliceAI delta; functional impact still observed in target tissue.

**Fix:** Use Pangolin per tissue or SpliceTransformer flags; require RNA validation in disease-relevant tissue.

### Pangolin: Masking and Out-of-Training Tissue

**Trigger:** `-m True`, or a disease tissue not represented in Pangolin's 4-species, 4-tissue (Cardoso-Moreira 2019 developmental) training set.

**Mechanism:** The mask erases real gains and losses (table above); Pangolin extrapolates poorly outside its training tissues.

**Symptom:** A splice-site loss or pseudoexon gain visible at `-m False` is 0.00 at `-m True`; score uncalibrated for the queried tissue and disagrees with patient RNA-seq.

**Fix:** Run `-m False`; fall back to SpliceAI for tissues not in the Pangolin training set, or run patient RNA-seq directly.

### MMSplice: Atypical Events

**Trigger:** Variant affecting a non-cassette event (MXE, complex multi-junction, AFE/ALE), or a deep-intronic variant.

**Mechanism:** MMSplice modular model is trained primarily on cassette exon events and scores variants only near an annotated exon.

**Symptom:** MMSplice ΔPSI doesn't match other predictors or empirical data, or the variant has no row.

**Fix:** Use SpliceAI for non-cassette events; restrict MMSplice to cassette exon contexts.

### CADD-Splice: Loss of Component Information

**Trigger:** Wanting to know which sub-component drove a high CADD-Splice score.

**Mechanism:** CADD-Splice combines splice predictors with the CADD annotations into a single C-score; sub-component contributions are abstracted.

**Symptom:** "High CADD-Splice score but unclear why."

**Fix:** Run SpliceAI and MMSplice separately to see which contributed.

### Branchpoint Variants: Low Information Motif

**Trigger:** Variant in the BPS region (-18 to -40 from 3'ss).

**Mechanism:** BPS motif (yUnAy) has low information content; CNNs struggle to learn the consensus.

**Symptom:** Confirmed BPS variant scores SpliceAI delta <0.2 despite functional disruption.

**Fix:** Require RNA validation (see Branchpoint Variant Detection).

## Population Database Lookup

| Database | Use for |
|----------|---------|
| gnomAD v4 | Allele frequency; SpliceAI annotations integrated |
| ClinVar | Existing classifications; SpliceAI integrated since 2020 |
| SpliceVarDB | Curated splice variants with experimental RNA validation |
| dbNSFP4 | Pre-computed splice scores aggregated |
| Recount3 | Tissue-specific PSI lookups from public RNA-seq |
| GTEx sQTL v8 | Tissue-specific splicing QTLs across 49 tissues |
| MaveDB | Splice MAVE results (e.g. BRCA1 saturation; Findlay 2018 *Nature*) |

Always check ClinVar first for existing classifications; cross-reference with gnomAD for population frequency before committing to PP3/PP4.

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ImportError: pkg_resources` / `No module named 'pkg_resources'` (SpliceAI, MMSplice, Pangolin) | setuptools >= 81 | `pip install "setuptools<81"` |
| SpliceAI/Pangolin skip a record: `Skipping record (ref issue)` / `Mismatch between FASTA (ref base: G) and variant file (ref base: A)`, exit 0 | REF does not match the FASTA (wrong build or coordinates) | Assert REF against the FASTA before scoring; report records missing from the output |
| Pangolin `Mismatch between FASTA (ref base: c) ...` on a correct REF | Soft-masked (lower-case) FASTA, e.g. UCSC `hg38.fa`; the comparison is case-sensitive | Use an upper-case FASTA (GENCODE genome) or `awk '/^>/{print;next}{print toupper($0)}'` |
| `ValueError: Duplicate ID CDS:ENST...` from `gffutils.create_db` | Real GENCODE GFF3 has repeated IDs | Use Pangolin's `create_db.py` on the GTF |
| `ModuleNotFoundError: pyfastx` / `vcf` / `torch` from Pangolin | `setup.py` declares no dependencies | Install them (see Install); `TypeError: Info.__new__() missing 'type_code'` -> `PyVCF3==1.0.0` |
| `pangolin: WARNING, skipping variant: Variant not contained in a gene body` | Contig names differ from the DB or gene absent from the (canonical) DB | Match chromosome naming; rebuild the DB from the full GTF |
| `Variant format not supported` (Pangolin) | REF and ALT both longer than one base with different lengths, or an allele with no ACGT (N-only) | Decompose into simple SNV/indel records; SpliceAI writes `.` scores for N-only alleles |
| SpliceAI `OSError: Can't write record` | Symbolic ALT such as `<DEL>` | Remove or expand structural variants |
| `mmsplice: Input sequence ... intron length cannot be longer than the input sequence` | FASTA truncated or not matching the GTF | Check FASTA size/index against the reference; use the GENCODE basic GTF |
| MMSplice output has no row for a variant | Variant not near an annotated exon (silent) | Compare output IDs with the input VCF |
| `SpliceVault` returns no rows | Variant not a catalogued splice-site SNV, or wrong REF/ALT/build | Check with SpliceAI which site is affected; GRCh38 only |
| `VariantValidator: invalid HGVS` | Wrong reference transcript or build | Specify NM_*.* version explicitly |
| A `chr`-prefixed VCF against a bare-contig FASTA (`chrX` vs `X`) | Naming mismatch | Both SpliceAI 1.3.1 and Pangolin 1.0.2 scored such records on the audit data; still keep names consistent (`bcftools annotate --rename-chrs`) if a tool skips |

## Common Pitfalls

- **Using SpliceAI score alone for clinical reporting** — combine with concordant predictors and ideally RNA validation; see ClinGen rules above.
- **Reading a missing score as benign** — `.`, no `SpliceAI=` tag, no Pangolin/MMSplice row = not scored.
- **Forgetting NMD direction** — confirmed splice disruption needs NMD-status check. Last-exon PTCs escape NMD and can be dominant-negative or gain-of-function.
- **Trusting LLMs for variant interpretation** — use as orchestrators on top of SpliceAI/VariantValidator/ClinVar; all clinical-grade calls require human expert sign-off.
- **Skipping HGVS validation** — invalid HGVS leads to silent reference-transcript mismatches; always run VariantValidator first.

## Related Skills

- splicing-qc - MaxEntScan + library QC for confirming predicted impact
- splicing-quantification - Empirical PSI from RNA-seq to validate predictions
- outlier-splicing-detection - FRASER2/DROP for RNA-seq confirmation in clinical samples
- variant-calling/clinical-interpretation - Broader ACMG/AMP variant interpretation framework
- variant-calling/variant-annotation - VEP plugin integration for SpliceAI

## References

- Jaganathan et al 2019 *Cell* - SpliceAI
- Zeng & Li 2022 *Genome Biol* - Pangolin
- Cheng et al 2019 *Genome Biol* - MMSplice
- Cheng et al 2021 *Genome Biol* - MTSplice (tissue MMSplice)
- You et al 2024 *Nat Commun* 15:9129 - SpliceTransformer
- Strauch et al 2022 *PLoS One* 17:e0269159 - CI-SpliceAI extended window
- Smith & Kitzman 2023 *Genome Biol* 24:294 - SpliceAI/Pangolin MPSA benchmark
- Rentzsch et al 2021 *Genome Med* - CADD-Splice
- Dawes et al 2023 *Nat Genet* - SpliceVault
- Walker et al 2023 *Am J Hum Genet* - ClinGen SVI splicing recommendations
- Riepe et al 2021 *Hum Mutat* 42:799 - SpliceAI in clinical pipelines (Riepe TV et al)
- Abou Tayoun et al 2018 *Hum Mutat* - PVS1 decision tree
- Richards et al 2015 *Genet Med* - ACMG/AMP framework
- den Dunnen et al 2016 *Hum Mutat* - HGVS standard
- Zhang et al 2022 *PNAS* (PMID 36306325) - BPHunter for branchpoints
- Hua et al 2008 *AJHG* - ISS-N1 / nusinersen mechanism
- Roberts et al 2023 *Nat Rev Drug Discov* 22:917-934 - DMD therapeutic approaches (exon-skipping ASOs)
- Findlay et al 2018 *Nature* - BRCA1 saturation genome editing (MAVE)
