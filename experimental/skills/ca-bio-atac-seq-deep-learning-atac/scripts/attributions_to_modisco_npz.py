#!/usr/bin/env python3
"""DeepLIFT/SHAP attributions from a chromBPNet accessibility model (chrombpnet_nobias.h5) -> modisco-lite inputs.

Loads the Keras .h5 in PyTorch (bpnet-lite BPNet.from_chrombpnet), runs tangermeme deep_lift_shap(hypothetical=True) on the counts or
profile head, and writes ohe.npz and attr.npz in the (N, 4, L) layout that `modisco motifs -s/-a` requires (central 1000 bp of each
2114 bp window). GPU is used when available.

usage: attributions_to_modisco_npz.py MODEL.h5 GENOME.fa PEAKS.narrowPeak OUTDIR [N_PEAKS=1500] [counts|profile]
  PEAKS: 10-column narrowPeak; windows are centred on start + summit (column 10). Peaks with N or off-chromosome windows are skipped.
env: dlatac-torch (torch, bpnet-lite, tangermeme, pyfaidx). Then:
  modisco motifs -s OUTDIR/ohe.npz -a OUTDIR/attr.npz -n 2000 -w 500 -o OUTDIR/modisco_results.h5
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
from bpnetlite.bpnet import BPNet, CountWrapper, ProfileWrapper
from pyfaidx import Fasta
from tangermeme.deep_lift_shap import deep_lift_shap
from tangermeme.utils import one_hot_encode

INPUT_LEN = 2114


def main(model_h5, genome, peaks, outdir, n_peaks=1500, head="counts"):
    device = "cuda" if torch.cuda.is_available() else "cpu"
    wrap = {"counts": CountWrapper, "profile": ProfileWrapper}[head]
    model = wrap(BPNet.from_chrombpnet(model_h5)).to(device).eval()
    fa = Fasta(genome)
    pk = pd.read_csv(peaks, sep="\t", header=None)
    if pk.shape[1] < 10:
        sys.exit("PEAKS must be 10-column narrowPeak (summit in column 10)")
    X = []
    for _, r in pk.head(int(n_peaks)).iterrows():
        s = int(r[1] + r[9]) - INPUT_LEN // 2
        if s < 0 or s + INPUT_LEN > len(fa[r[0]]):
            continue
        seq = str(fa[r[0]][s:s + INPUT_LEN]).upper()
        if set(seq) <= set("ACGT"):
            X.append(one_hot_encode(seq))
    X = torch.stack(X)
    A = deep_lift_shap(model, X, hypothetical=True, n_shuffles=10, batch_size=32, device=device, random_state=0, verbose=False)
    lo = (INPUT_LEN - 1000) // 2
    os.makedirs(outdir, exist_ok=True)
    np.savez_compressed(os.path.join(outdir, "ohe.npz"), X[:, :, lo:lo + 1000].numpy().astype(np.int8))
    np.savez_compressed(os.path.join(outdir, "attr.npz"), A[:, :, lo:lo + 1000].numpy().astype(np.float32))
    print("wrote %d peaks x %s head to %s (ohe.npz, attr.npz; shape (N,4,1000))" % (len(X), head, outdir))


if __name__ == "__main__":
    if not 5 <= len(sys.argv) <= 7:
        sys.exit(__doc__)
    main(*sys.argv[1:])
