#!/usr/bin/env python
"""SYNTHETIC cohort for auditing bio-outlier-splicing-detection.  Nothing here is real data.

Simulates paired-end (2x75) RNA-seq alignments for N_SAMPLES samples over N_GENES multi-exon genes on a fake
contig `chr21`, writes coordinate-sorted, indexed BAMs (STAR-like tags NH/HI/XS), a truth table of PLANTED
aberrant-splicing / expression events, and a gene x sample fragment-count matrix.

Planted events (truth in planted_truth.tsv):
  S05  gene g010  exon-2 skipping 0.05 -> 0.75           (annotated-junction usage outlier)
  S12  gene g030  cryptic donor, exon1 extended +45 nt   (novel 5'ss, annotation-free)
  S20  gene g050  intron-1 retention 0.03 -> 0.70        (intron retention, unspliced reads)
  S25  gene g070  100-nt pseudoexon in intron 1, 0.5     (two novel junctions)
  S15  gene g090  expression x0.05 (LoF/NMD-like)        (OUTRIDER-style, splicing unchanged)
  S28  gene g100  expression x6                          (over-expression)
  S29  tissue-mismatch sample: half of all genes have 5x higher baseline skipping (global shift; no planted gene)
Run:  python 01_make_synth_cohort.py <outdir>  (uses pysam; run in WSL env as-core)
"""
import sys, os, json
import numpy as np
import pysam

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
N_SAMPLES, N_GENES = 30, 200
CONTIG, CLEN = "chr21", 2_000_000
rng = np.random.default_rng(20260920)
RL = 75

# ---- gene models -------------------------------------------------------------------------------
genes = []
for i in range(N_GENES):
    ne = int(rng.integers(3, 6))
    start = 10_000 + i * 8_000
    ex, pos = [], start
    for k in range(ne):
        L = int(rng.integers(150, 251))
        ex.append((pos, pos + L))
        pos += L + int(rng.integers(600, 1500))
    genes.append(dict(id=f"g{i:03d}", ex=ex, mu=float(np.exp(rng.normal(np.log(320), 0.5))),
                      base_skip=float(rng.beta(2, 30)) if rng.random() < 0.5 else 0.0,
                      base_ir=0.02))
gid = {g["id"]: g for g in genes}
sid = [f"S{j+1:02d}" for j in range(N_SAMPLES)]
depth = np.exp(rng.normal(0, 0.25, N_SAMPLES))
batch = (np.arange(N_SAMPLES) % 2)             # two mild batches: batch 1 has 1.5x baseline skipping
mismatch_sample, mismatch_genes = "S29", set(g["id"] for g in genes[::2])

planted = [
    dict(sample="S05", gene="g010", type="exon_skipping"),
    dict(sample="S12", gene="g030", type="cryptic_donor"),
    dict(sample="S20", gene="g050", type="intron_retention"),
    dict(sample="S25", gene="g070", type="pseudoexon"),
    dict(sample="S15", gene="g090", type="expression_down"),
    dict(sample="S28", gene="g100", type="expression_up"),
]
for p in planted:
    g = gid[p["gene"]]; p["chrom"] = CONTIG; p["gene_start"] = g["ex"][0][0]; p["gene_end"] = g["ex"][-1][1]
with open(os.path.join(out, "planted_truth.tsv"), "w") as f:
    f.write("sample\tgene\ttype\tchrom\tgene_start\tgene_end\n")
    for p in planted:
        f.write(f"{p['sample']}\t{p['gene']}\t{p['type']}\t{p['chrom']}\t{p['gene_start']}\t{p['gene_end']}\n")
with open(os.path.join(out, "genes.tsv"), "w") as f:
    f.write("gene\tchrom\tstart\tend\tstrand\texons\n")
    for g in genes:
        f.write(f"{g['id']}\t{CONTIG}\t{g['ex'][0][0]}\t{g['ex'][-1][1]}\t+\t{','.join(f'{a}-{b}' for a,b in g['ex'])}\n")

def merge(segs):
    out_ = [list(segs[0])]
    for a, b in segs[1:]:
        if a == out_[-1][1]: out_[-1][1] = b
        else: out_.append([a, b])
    return [tuple(s) for s in out_]

def isoform_segments(g, kind, extra=None):
    ex = list(g["ex"])
    if kind == "canonical":
        segs = ex
    elif kind == "skip":                       # skip exon index 1
        segs = [ex[0]] + ex[2:]
    elif kind == "retain":                     # retain intron `extra`
        k = extra; segs = ex[:k] + [(ex[k][0], ex[k + 1][1])] + ex[k + 2:]
        segs = [tuple(s) for s in segs]
    elif kind == "cryptic":                    # exon 1 donor +45 nt
        segs = [(ex[0][0], ex[0][1] + 45)] + ex[1:]
    elif kind == "pseudo":                     # pseudoexon (100 nt) in the middle of intron 1
        mid = (ex[0][1] + ex[1][0]) // 2
        segs = [ex[0], (mid, mid + 100)] + ex[1:]
    return merge(segs)

def map_read(segs, cum, off, ln):
    """genomic pos (0-based) + cigar of transcript interval [off, off+ln)"""
    ops, pos0, remaining, cur = [], None, ln, off
    for (a, b), c0 in zip(segs, cum):
        c1 = c0 + (b - a)
        if cur >= c1: continue
        take = min(remaining, c1 - cur)
        gs = a + (cur - c0)
        if pos0 is None: pos0 = gs
        else:
            gap = gs - last_end
            ops.append((3, gap))
        ops.append((0, take)); last_end = gs + take
        cur += take; remaining -= take
        if remaining == 0: break
    return pos0, ops

header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": CONTIG, "LN": CLEN}]}
counts = np.zeros((N_GENES, N_SAMPLES), dtype=int)
qual = pysam.qualitystring_to_array("I" * RL)

for j, s in enumerate(sid):
    recs = []
    for gi, g in enumerate(genes):
        mu = g["mu"] * depth[j]
        if s == "S15" and g["id"] == "g090": mu *= 0.05
        if s == "S28" and g["id"] == "g100": mu *= 6.0
        n = int(rng.negative_binomial(20, 20 / (20 + mu)))
        counts[gi, j] = n
        # isoform probabilities for this gene x sample
        bs = g["base_skip"] * (1.5 if batch[j] else 1.0)
        if s == mismatch_sample and g["id"] in mismatch_genes: bs = max(bs, 0.02) * 5
        bs = min(0.9, float(rng.beta(max(bs, 1e-3) * 80 + 1e-3, (1 - min(bs, 0.95)) * 80))) if bs > 0 else 0.0
        ir = g["base_ir"]; pcrypt = 0.0; ppseudo = 0.0
        if s == "S05" and g["id"] == "g010": bs = 0.75
        if s == "S20" and g["id"] == "g050": ir = 0.70
        if s == "S12" and g["id"] == "g030": pcrypt = 0.6
        if s == "S25" and g["id"] == "g070": ppseudo = 0.5
        ne = len(g["ex"])
        iso_cache = {}
        for _ in range(n):
            u = rng.random()
            if u < pcrypt: key = ("cryptic", None)
            elif u < pcrypt + ppseudo: key = ("pseudo", None)
            elif u < pcrypt + ppseudo + ir: key = ("retain", 0 if (s == "S20" and g["id"] == "g050") else int(rng.integers(0, ne - 1)))
            elif u < pcrypt + ppseudo + ir + bs: key = ("skip", None)
            else: key = ("canonical", None)
            if key not in iso_cache:
                segs = isoform_segments(g, key[0], key[1])
                cum = np.concatenate([[0], np.cumsum([b - a for a, b in segs])])[:-1]
                iso_cache[key] = (segs, cum, sum(b - a for a, b in segs))
            segs, cum, tl = iso_cache[key]
            L = int(min(tl, rng.integers(200, 321)))
            if L < 2 * RL - 20: continue
            f0 = int(rng.integers(0, tl - L + 1))
            p1, c1 = map_read(segs, cum, f0, RL)
            p2, c2 = map_read(segs, cum, f0 + L - RL, RL)
            recs.append((p1, c1, p2, c2, len(recs)))
    recs.sort(key=lambda r: r[0])
    path = os.path.join(out, f"{s}.unsorted.bam")
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for (p1, c1, p2, c2, k) in recs:
            name = f"{s}_r{k}"
            tlen = (p2 + sum(l for o, l in c2 if o in (0, 3))) - p1
            for first, (p, c) in ((True, (p1, c1)), (False, (p2, c2))):
                a = pysam.AlignedSegment(bam.header)
                a.query_name = name; a.query_sequence = "A" * RL; a.query_qualities = qual
                a.flag = 99 if first else 147
                a.reference_id = 0; a.reference_start = p; a.mapping_quality = 255
                a.cigartuples = c
                a.next_reference_id = 0; a.next_reference_start = p2 if first else p1
                a.template_length = tlen if first else -tlen
                a.set_tags([("NH", 1, "i"), ("HI", 1, "i"), ("XS", "+", "A")])
                bam.write(a)
    final = os.path.join(out, f"{s}.bam")
    pysam.sort("-o", final, path); pysam.index(final); os.remove(path)
    print(s, "fragments", len(recs), flush=True)

with open(os.path.join(out, "gene_counts.tsv"), "w") as f:
    f.write("gene\t" + "\t".join(sid) + "\n")
    for gi, g in enumerate(genes):
        f.write(g["id"] + "\t" + "\t".join(map(str, counts[gi])) + "\n")
with open(os.path.join(out, "README_SYNTHETIC.txt"), "w") as f:
    f.write("SYNTHETIC data generated by 01_make_synth_cohort.py (seed 20260920). Not real. Truth: planted_truth.tsv\n")
print("done")
