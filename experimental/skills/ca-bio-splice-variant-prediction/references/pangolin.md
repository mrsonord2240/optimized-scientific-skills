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

**Tissue-specific output.** The CLI prints the maximum gain and maximum loss over the four tissue models (heart, liver, brain, testis), not per-tissue scores. Per-tissue values come from `scripts/pangolin_tissue.py`, which reuses the installed models:

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
