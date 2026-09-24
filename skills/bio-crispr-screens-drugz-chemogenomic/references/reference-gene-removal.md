# Removing Reference Genes from the Analysis

**Goal:** Keep reference essential or control genes from inflating the Z-score null distribution.

**Approach:** `-r` takes a **comma-delimited list of gene symbols, not a file path** (`drugz.py` does
`args.remove_genes.split(',')`). Passing a filename removes nothing and still exits 0. The named genes
are dropped from the input before any Z-scoring, so they also disappear from the output file: `-r` is
an exclusion, not a re-weighting.

```bash
# Inline list
python drugz.py \
    -i counts.txt \
    -o drugz_clean.txt \
    -c Veh_r1,Veh_r2 \
    -x Drug_r1,Drug_r2 \
    -r RPS3,RPL11,EIF3A,POLR2A,CDK1
```

```bash
# From a reference set: take the first column, skip the header, join with commas.
# CEGv2.txt (hart-lab/bagel) is tab-separated with a header (GENE, HGNC_ID, ENTREZ_ID);
# joining its raw lines gives tokens like "AARS<TAB>HGNC:20<TAB>16", which match no gene
# and silently exclude nothing.
CEG=$(tail -n +2 CEGv2.txt | cut -f1 | paste -sd, -)
python drugz.py -i counts.txt -o drugz_clean.txt -c Veh_r1,Veh_r2 -x Drug_r1,Drug_r2 -r "$CEG"
```

**Verify the exclusion happened** -- the tool cannot tell you it matched nothing:

```bash
# genes in the unfiltered output but not the filtered one; should equal the number you excluded
comm -23 <(cut -f1 drugz_output.txt | sort) <(cut -f1 drugz_clean.txt | sort) | wc -l
```

**When to use:** If pilot drugZ runs show many essential genes appearing as "sensitizers" purely
because they drop out under any condition, removing them gives a cleaner drug-specific signal.
