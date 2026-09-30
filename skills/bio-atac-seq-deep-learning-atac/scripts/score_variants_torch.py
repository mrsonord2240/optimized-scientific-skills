#!/usr/bin/env python3
"""Score SNPs with a chromBPNet accessibility model (chrombpnet_nobias.h5) in PyTorch on GPU or CPU.

Loads the Keras .h5 through bpnet-lite (BPNet.from_chrombpnet, no TensorFlow), uses tangermeme substitution_effect on the
LOG-COUNT head and reports log2FC = (log_count_alt - log_count_ref) / ln 2, the quantity variant-scorer reports as `logfc`.
Forward strand only: it matches variant-scorer run with --forward_only (its default also averages the reverse complement).

usage: score_variants_torch.py MODEL.h5 GENOME.fa VARIANTS.tsv OUT.tsv
  VARIANTS.tsv: 5 tab-separated headerless columns  chr pos(1-based) ref alt variant_id   (SNVs only)
env: dlatac-torch (torch, bpnet-lite, tangermeme, pyfaidx)
"""
import sys

import numpy as np
import pandas as pd
import torch
from bpnetlite.bpnet import BPNet, CountWrapper
from pyfaidx import Fasta
from tangermeme.utils import one_hot_encode
from tangermeme.variant_effect import substitution_effect

BASES = "ACGT"


def main(model_h5, genome, variants, out):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = CountWrapper(BPNet.from_chrombpnet(model_h5)).to(device).eval()   # CountWrapper: log-count head only
    fa = Fasta(genome)
    v = pd.read_csv(variants, sep="\t", header=None, names=["chr", "pos", "ref", "alt", "id"])
    half = 1057  # ChromBPNet input length 2114; variant sits at index 1057
    keep, seqs = [], []
    for i, r in v.iterrows():
        s = int(r.pos) - 1 - half
        if s < 0 or s + 2 * half > len(fa[r.chr]):
            continue
        seq = str(fa[r.chr][s:s + 2 * half]).upper()
        if seq[half] != r.ref or r.alt not in BASES or "N" in seq:
            continue  # reference mismatch, non-ACGT allele or gap: skipped, counted below
        keep.append(i)
        seqs.append(one_hot_encode(seq))
    if not keep:
        sys.exit("no scorable variants (check chromosome names, 1-based positions and REF alleles)")
    X = torch.stack(seqs)                                                     # (N, 4, 2114)
    subs = torch.tensor([[j, half, BASES.index(v.alt[i])] for j, i in enumerate(keep)])
    y_ref, y_alt = substitution_effect(model, X, subs, device=device, verbose=False)
    y_ref, y_alt = np.asarray(y_ref).ravel(), np.asarray(y_alt).ravel()
    res = v.loc[keep].copy()
    res["logcount_ref"], res["logcount_alt"] = y_ref, y_alt
    res["log2fc"] = (y_alt - y_ref) / np.log(2)
    res.to_csv(out, sep="\t", index=False)
    print("scored %d of %d variants (%d skipped) -> %s" % (len(keep), len(v), len(v) - len(keep), out))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
