#!/usr/bin/env python
"""Intersect thresholded ABC and ENCODE-rE2G links and optionally flag HiChIP loop support.

Both inputs are the pipelines' own thresholded tables, so each method's calibrated cut-off is
already applied. Links are matched on (`name`, `TargetGene`): `name` encodes the element
coordinates, so this is valid only when both methods used the same candidate regions (run
ENCODE-rE2G on the ABC peak set, or match by interval overlap instead).

HiChIP support: a link is supported when a loop has one anchor overlapping the enhancer and the
other overlapping the gene TSS (either orientation). Loops are a BEDPE (chr1 s1 e1 chr2 s2 e2 ...).

Usage: combine_predictions.py --abc EnhancerPredictionsFull_threshold*.tsv \
           --re2g encode_e2g_predictions_threshold*.tsv.gz [--loops loops.bedpe] --out high_confidence.tsv
"""
import argparse

import pandas as pd
import pyranges as pr


def loop_support(links, loops):
    """Return a boolean Series (indexed like links) marking links spanned by a loop."""
    enh = pr.PyRanges(pd.DataFrame({"Chromosome": links["chr"], "Start": links["start"],
                                    "End": links["end"], "link": links.index}))
    tss = pr.PyRanges(pd.DataFrame({"Chromosome": links["chr"], "Start": links["TargetGeneTSS"].astype(int),
                                    "End": links["TargetGeneTSS"].astype(int) + 1, "link": links.index}))
    hits = set()
    for a, b in ((("chr1", "s1", "e1"), ("chr2", "s2", "e2")), (("chr2", "s2", "e2"), ("chr1", "s1", "e1"))):
        def anchors(cols):
            return pr.PyRanges(pd.DataFrame({"Chromosome": loops[cols[0]], "Start": loops[cols[1]],
                                             "End": loops[cols[2]], "loop": loops.index}))
        e = enh.join(anchors(a)).df[["link", "loop"]]
        t = tss.join(anchors(b)).df[["link", "loop"]]
        hits |= set(map(tuple, e.merge(t, on=["link", "loop"])[["link", "loop"]].values))
    return links.index.isin({link for link, _ in hits})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--abc", required=True)
    ap.add_argument("--re2g", required=True)
    ap.add_argument("--loops")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    abc = pd.read_csv(args.abc, sep="\t")
    re2g = pd.read_csv(args.re2g, sep="\t")
    key = ["name", "TargetGene"]
    both = abc.merge(re2g[key + ["ENCODE-rE2G.Score"]], on=key)
    print(f"ABC {len(abc)} links, ENCODE-rE2G {len(re2g)} links, both {len(both)}")
    if args.loops:
        loops = pd.read_csv(args.loops, sep="\t", header=None, usecols=range(6),
                            names=["chr1", "s1", "e1", "chr2", "s2", "e2"], comment="#")
        both["hichip_support"] = loop_support(both, loops)
        print(f"{int(both['hichip_support'].sum())} of {len(both)} supported by HiChIP loops")
    both.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
