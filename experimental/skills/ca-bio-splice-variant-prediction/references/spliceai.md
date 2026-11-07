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

`-D` = **maximum distance between the variant and the gained/lost splice site** (default 50). `-M 0` (default) returns raw scores; `-M 1` masks splice gains at annotated sites and losses at unannotated sites (annotation = SpliceAI's canonical GENCODE table). Output INFO format: `SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL`. Delta score = max(DS_AG, DS_AL, DS_DG, DS_DL). DS labels: AG = acceptor gain, AL = acceptor loss, DG = donor gain, DL = donor loss; DP = position of that site relative to the variant. One record can carry several annotations (readthrough/overlapping genes such as RPL36A-HNRNPH2 add rows): `spliceai_clingen_classify.py` reports both `top_score_gene` and `annotated_genes`; neither selects the disease-relevant transcript. Resolve the MANE/clinical transcript independently.

**SpliceAI exits 0 on records it cannot score.** REF mismatch, a deletion longer than 2*D, or no overlapping gene give a record with **no `SpliceAI=` tag** (stderr warning only); N-containing sequence gives `SpliceAI=...|.|.|.|.|.|.|.|.` (all dots). A `<DEL>` symbolic ALT crashes it. The shipped classifier pre-filters `*` and symbolic ALT alleles so they appear as `not_scored` warnings while supported alleles in the same batch still run. Assert REF against the FASTA first and compare output to input: a failed variant must not read as benign or vanish.

Parse, label and check with `scripts/splice_parsers.py` (`.` -> NaN -> `not_scored`, boundaries inclusive, unmatched records reported):

```python
from scripts.splice_parsers import parse_spliceai_vcf, classify_delta, read_input_vcf, unscored_report

s = parse_spliceai_vcf('output.vcf')                       # one row per variant x ALT x gene
s['acmg_evidence'] = classify_delta(s['delta_max'])        # BP4 / inconclusive / PP3_supporting[_prec0.5|_prec0.8] / not_scored
print(unscored_report('input.vcf', s['key'], 'SpliceAI'))  # input records with no SpliceAI= tag
```

`scripts/spliceai_clingen_classify.py input.vcf genome.fa --build grch38` runs the whole flow (default window, wider re-run of everything scoring < 0.20 including 0.00, PP3/BP4 labels, skipped-record report). `python test_splice_parsers.py` in `examples/` self-tests the parsers and boundary labels on real SpliceAI/Pangolin/MMSplice output (`test_data/`, GRCh38 panel: TP53 c.673-2A>G, DMD c.9563+1G>A / c.31+1G>A, OTC c.386+5G>A, GLA c.370-1G>A, GLA c.639+919G>A, three ClinVar-benign GLA variants). Expected SpliceAI D50 delta_max on that panel: TP53 1.00, DMD c.9563+1G>A 0.99, DMD c.31+1G>A 0.95, OTC c.386+5G>A 0.94, GLA c.370-1G>A 1.00, GLA c.639+919G>A 0.30, the three benign GLA variants 0.00.
