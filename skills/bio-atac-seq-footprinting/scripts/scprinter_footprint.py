#!/usr/bin/env python3
"""Multi-scale footprint scores with scPrinter 1.2.0 (classic get_footprint_score route), bulk or per cell group.

Steps: build a custom Genome, predict the genome-wide Tn5 bias (GPU strongly advised), import one fragment file, then score
every region in --regions at each footprint scale (mode), for one pseudobulk (default) or for each group in --groups.

Usage:
  scprinter_footprint.py --fragments frags.tsv.gz --fasta hg38.fa --gtf genes.gtf --blacklist bl.bed \
      --regions sites.bed --outdir out [--groups cells.tsv] [--shift 0,0|4,-5|auto] [--bias bias.h5] [--modes 2-100] \
      [--width 200] [--device cuda:0]

--fragments  bgzip + tabix BED with columns chrom, start, end, barcode.
--shift      Tn5 shift ALREADY applied to the fragment ends (scPrinter's plus_shift, minus_shift): 0,0 (default) for raw
             fragments from a BAM (usage-guide recipe), 4,-5 for Cell Ranger fragments; 'auto' uses scPrinter's beta detection.
--groups     TSV without header: barcode, group (for example a cluster label). Only listed barcodes are imported and each
             group is scored as one pseudobulk. Without it the whole file is one pseudobulk sample.
--regions    BED; each region is resized to --width around its centre (motif centres for footprint work).
Outputs: <outdir>/footprints.npz (scores: regions x modes x width, or regions x groups x modes x width with --groups;
keys, groups, modes) and <outdir>/center_by_mode.tsv (group, mode, centre and flank means).
Exit 2 when an input file is missing. Set SCPRINTER_DATA to a writable directory: the first import downloads the pretrained
models (needs internet).
"""
import argparse
import os
import sys


def parse_modes(text):
    import numpy as np
    if "-" in text:
        lo, hi = (int(x) for x in text.split("-"))
        return np.arange(lo, hi + 1)
    return np.array([int(x) for x in text.split(",")])


def parse_shift(text):
    if text == "auto":
        return True, 4, -5
    try:
        plus, minus = (int(x) for x in text.split(","))
    except ValueError:
        print(f"ERROR: --shift must be 'auto' or 'PLUS,MINUS', got {text!r}", file=sys.stderr)
        sys.exit(2)
    return False, plus, minus


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fragments", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--blacklist", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--groups", help="TSV barcode<TAB>group; score each group as a pseudobulk")
    ap.add_argument("--shift", default="0,0", help="shift already applied to fragment ends: 'PLUS,MINUS' or 'auto' (default 0,0)")
    ap.add_argument("--bias", help="existing genome Tn5 bias .h5; predicted into <outdir>/bias.h5 when absent")
    ap.add_argument("--modes", default="2-100", help="scales: 'lo-hi' or comma list (default 2-100)")
    ap.add_argument("--width", type=int, default=200)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--name", default="sample")
    ap.add_argument("--min-fragments", type=int, default=1000, help="per-barcode fragment minimum at import")
    a = ap.parse_args()

    inputs = {"--fragments": a.fragments, "--fasta": a.fasta, "--gtf": a.gtf, "--blacklist": a.blacklist,
              "--regions": a.regions, "--bias": a.bias, "--groups": a.groups}
    for flag, path in inputs.items():
        if path is not None and not os.path.isfile(path):
            print(f"ERROR: input not found: {flag} {path}", file=sys.stderr)
            sys.exit(2)
    auto_shift, plus_shift, minus_shift = parse_shift(a.shift)

    os.makedirs(a.outdir, exist_ok=True)
    os.environ.setdefault("SCPRINTER_DATA", os.path.join(os.path.abspath(a.outdir), "scprinter_data"))
    os.makedirs(os.environ["SCPRINTER_DATA"], exist_ok=True)
    import numpy as np
    import pandas as pd
    import scprinter as scp

    modes = parse_modes(a.modes)
    bias = a.bias or os.path.join(a.outdir, "bias.h5")
    if not os.path.exists(bias):
        scp.genome.predict_genome_tn5_bias(fa_file=a.fasta, save_name=bias, tn5_model=scp.datasets.pretrained_Tn5_bias_model,
                                           context_radius=50, device=a.device, batch_size=5000)
    genome = scp.genome.Genome(name="custom", fa_file=a.fasta, gff_file=a.gtf, bias_file=bias, blacklist_file=a.blacklist)
    h5 = os.path.join(a.outdir, "printer.h5ad")
    if os.path.exists(h5):
        os.remove(h5)
    shift = dict(auto_detect_shift=auto_shift, plus_shift=plus_shift, minus_shift=minus_shift)
    if a.groups:
        cells = pd.read_csv(a.groups, sep="\t", header=None, dtype=str).iloc[:, :2]
        cells.columns = ["barcode", "group"]
        printer = scp.pp.import_fragments(path_to_frags=a.fragments, barcodes=list(cells.barcode), savename=h5, genome=genome,
                                          min_num_fragments=a.min_fragments, min_tsse=0, sorted_by_barcode=False,
                                          low_memory=False, **shift)
        kept = set(printer.obs_names)
        names = sorted(cells.group.unique())
        groups = [[b for b in cells.barcode[cells.group == g] if b in kept] for g in names]
        print("cells per group after import:", dict(zip(names, map(len, groups))), f"({len(cells) - len(kept)} dropped)")
        if any(len(g) == 0 for g in groups):
            print("ERROR: a group has no imported cells; lower --min-fragments or check the barcodes", file=sys.stderr)
            sys.exit(1)
    else:
        printer = scp.pp.import_fragments(path_to_frags=[a.fragments], barcodes=[None], savename=h5, genome=genome,
                                          sample_names=[a.name], min_num_fragments=a.min_fragments, min_tsse=0,
                                          sorted_by_barcode=False, low_memory=False, **shift)
        names, groups = [a.name], [list(printer.obs_names)]
    printer.load_disp_model()
    regions = pd.read_csv(a.regions, sep="\t", header=None, comment="#").iloc[:, :3]
    regions.columns = ["chrom", "start", "end"]
    scp.tl.get_footprint_score(printer, groups, names, regions, region_width=a.width, modes=modes,
                               footprintRadius=None, flankRadius=None, n_jobs=a.jobs, save_key="fp", backed=True,
                               overwrite=True)
    ad = printer.footprintsadata["fp"]
    keys = list(ad.obsm.keys())                           # "chrom:start-end" of the resized regions, same order as scores
    scores = np.stack([np.asarray(ad.obsm[k]) for k in keys])      # regions x groups x modes x width
    printer.close()
    if not np.isfinite(scores).all():
        print("ERROR: non-finite footprint scores", file=sys.stderr)
        sys.exit(1)
    mid = a.width // 2
    rows = []
    for gi, g in enumerate(names):
        s = scores[:, gi]
        center = s[:, :, mid - 10:mid + 10].mean(axis=2).mean(axis=0)
        flank = np.concatenate([s[:, :, :mid - 40], s[:, :, mid + 40:]], axis=2).mean(axis=2).mean(axis=0)
        rows.append(pd.DataFrame({"group": g, "mode": modes, "center_pm10bp_mean": center, "flank_mean": flank}))
    pd.concat(rows).to_csv(os.path.join(a.outdir, "center_by_mode.tsv"), sep="\t", index=False)
    if not a.groups:
        scores = scores[:, 0]
    np.savez_compressed(os.path.join(a.outdir, "footprints.npz"), scores=scores, keys=np.array(keys),
                        groups=np.array(names), modes=modes)
    print(f"scprinter {scp.__version__}; scores {scores.shape}; wrote {a.outdir}/footprints.npz")


if __name__ == "__main__":
    main()
