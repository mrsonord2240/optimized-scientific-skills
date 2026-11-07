#!/usr/bin/env python3
"""Enformer (enformer-pytorch, weights EleutherAI/enformer-official-rough from Hugging Face) SNP effect on chosen output tracks.

Each SNP gets a 196,608 bp window centred on it (Enformer's fixed input length; output = 896 bins x 128 bp x 5313 human tracks,
covering the central 114,688 bp). Effect = log2((alt + 1) / (ref + 1)) of the summed prediction in the two 128 bp bins that flank
the variant. Forward strand only, no shift/reverse-complement averaging, pseudocount 1 in prediction units.

usage: enformer_variant_effect.py GENOME.fa VARIANTS.tsv TRACKS OUT.tsv
  VARIANTS.tsv  5 tab-separated headerless columns: chr pos(1-based) ref alt variant_id (SNVs)
  TRACKS        comma-separated 0-based track indices from targets_human.txt
                (https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt, column `index`;
                 e.g. 12 = DNASE:GM12878). The model file carries no track names.
env: dlatac-torch (torch, enformer-pytorch, pyfaidx); ~2 GB weight download on first use, GPU strongly preferred.
"""
import sys

import numpy as np
import pandas as pd
import torch
from enformer_pytorch import Enformer
from pyfaidx import Fasta

LEN = 196_608
CODE = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
CENTER_BINS = (447, 449)  # output bins 447 and 448 flank the window centre (bin 448 starts at the variant)


def encode(seq):
    return torch.tensor([CODE[c] for c in seq])


def main(genome, variants, tracks, out):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    tracks = [int(t) for t in tracks.split(",")]
    model = Enformer.from_pretrained("EleutherAI/enformer-official-rough").to(device).eval()
    fa = Fasta(genome)
    v = pd.read_csv(variants, sep="\t", header=None, names=["chr", "pos", "ref", "alt", "id"])
    rows = []
    for _, r in v.iterrows():
        s = int(r.pos) - 1 - LEN // 2
        if s < 0 or s + LEN > len(fa[r.chr]):
            continue
        ref = str(fa[r.chr][s:s + LEN]).upper()
        if ref[LEN // 2] != r.ref or r.alt not in "ACGT":
            continue
        alt = ref[:LEN // 2] + r.alt + ref[LEN // 2 + 1:]
        with torch.no_grad():
            y = model(torch.stack([encode(ref), encode(alt)]).to(device))["human"]      # (2, 896, 5313)
        c = y[:, CENTER_BINS[0]:CENTER_BINS[1], tracks].sum(1).cpu().numpy()           # (2, n_tracks)
        rows.append([r.chr, r.pos, r.ref, r.alt, r.id] + list(np.log2((c[1] + 1) / (c[0] + 1))))
    res = pd.DataFrame(rows, columns=["chr", "pos", "ref", "alt", "id"] + ["log2fc_track%d" % t for t in tracks])
    res.to_csv(out, sep="\t", index=False)
    print("scored %d of %d variants -> %s" % (len(res), len(v), out))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
