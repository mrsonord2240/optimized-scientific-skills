#!/usr/bin/env python3
'''
Per-tissue Pangolin scores for one variant. The `pangolin` CLI reports only the maximum over the four tissue
models; this reuses the installed package's models and encoding to print heart, liver, brain and testis
separately (same arithmetic as pangolin.pangolin.compute_score, minus the max over tissues).

Usage: python pangolin_tissue.py chr17 7674292 T C --strand - --fasta genome.fa [-d 50]
`--strand` is the strand of the gene (Pangolin takes it from the annotation DB; here you supply it).
No annotation masking is applied (equivalent to `-m False`). Substitutions and simple indels only.
Research-use decision support. Checked on Pangolin 1.0.2 (torch 2.13 CPU); the tissue maximum reproduces
the CLI's `-m False` output (TP53 c.673-2A>G: +0.72 at +47, -0.89 at -2).
'''
import argparse
import numpy as np
import pyfastx
import torch
from pkg_resources import resource_filename
from pangolin.model import L, W, AR, Pangolin
from pangolin.pangolin import one_hot_encode

TISSUES = ['heart', 'liver', 'brain', 'testis']   # model numbers 0, 2, 4, 6; 3 replicate models each; channels 1, 4, 7, 10


def load_models():
    models = []
    for i in [0, 2, 4, 6]:
        for j in range(1, 4):
            m = Pangolin(L, W, AR)
            m.load_state_dict(torch.load(resource_filename('pangolin', 'models/final.%s.%s.3.v2' % (j, i)),
                                         map_location=torch.device('cpu')))
            m.eval()
            models.append(m)
    return models


def tissue_scores(ref_seq, alt_seq, strand, d, models):
    '''(4, 2d+1) array of alt-minus-ref splice-site score per tissue, position -d..+d around the variant.'''
    ref_t = torch.from_numpy(np.expand_dims(one_hot_encode(ref_seq, strand).T, 0)).float()
    alt_t = torch.from_numpy(np.expand_dims(one_hot_encode(alt_seq, strand).T, 0)).float()
    out = []
    for j in range(4):
        per_model = []
        for model in models[3 * j:3 * j + 3]:
            with torch.no_grad():
                ref = model(ref_t)[0][[1, 4, 7, 10][j], :].numpy()
                alt = model(alt_t)[0][[1, 4, 7, 10][j], :].numpy()
            if strand == '-':
                ref, alt = ref[::-1], alt[::-1]
            l, nd = 2 * d + 1, abs(len(ref) - len(alt))
            if len(ref) > len(alt):
                alt = np.concatenate([alt[0:l // 2 + 1], np.zeros(nd), alt[l // 2 + 1:]])
            elif len(ref) < len(alt):
                alt = np.concatenate([alt[0:l // 2], np.max(alt[l // 2:l // 2 + nd + 1], keepdims=True), alt[l // 2 + nd + 1:]])
            per_model.append(alt - ref)
        out.append(np.mean(per_model, axis=0))
    return np.array(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('chrom'); ap.add_argument('pos', type=int); ap.add_argument('ref'); ap.add_argument('alt')
    ap.add_argument('--strand', required=True, choices=['+', '-']); ap.add_argument('--fasta', required=True)
    ap.add_argument('-d', '--distance', type=int, default=50)
    a = ap.parse_args()
    d = a.distance
    fa = pyfastx.Fasta(a.fasta)
    chrom = a.chrom if a.chrom in fa.keys() else ('chr' + a.chrom if 'chr' + a.chrom in fa.keys() else a.chrom[3:])
    seq = fa[chrom][a.pos - 5001 - d:a.pos + len(a.ref) + 4999 + d].seq
    if seq[5000 + d:5000 + d + len(a.ref)] != a.ref:
        raise SystemExit('REF mismatch: FASTA has %s at %s:%d (soft-masked FASTA? use an upper-case one)'
                         % (seq[5000 + d:5000 + d + len(a.ref)], chrom, a.pos))
    alt_seq = seq[:5000 + d] + a.alt + seq[5000 + d + len(a.ref):]
    sc = tissue_scores(seq, alt_seq, a.strand, d, load_models())
    for t, s in zip(TISSUES, sc):
        print(f'{t}: gain {s.max():+.2f} at {int(s.argmax()) - d:+d}, loss {s.min():+.2f} at {int(s.argmin()) - d:+d}')
    print(f'tissue maximum (what the CLI reports): gain {sc.max():+.2f}, loss {sc.min():+.2f}')
