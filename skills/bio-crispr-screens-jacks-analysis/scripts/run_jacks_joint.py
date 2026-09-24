# Purpose: run JACKS jointly across screens through the Python API (jacks.jacks_io.runJACKS).
# Inputs:  count matrix, replicate map (Replicate, Sample, Control; tab-separated WITH header) and
#          guide map (sgRNA, Gene; tab-separated WITH header). Checked on JACKS 0.2.
# Usage:   python run_jacks_joint.py counts.txt replicatemap.txt guidemap.txt --outprefix jacks_out
#              [--ctrl-genes NEGv1.txt --n-pseudo 2000 --seed 1] [--apply-w-hp]
#          For reproducible p-values run as: PYTHONHASHSEED=1 python run_jacks_joint.py ... --seed 1
#          Use --common-ctrl-sample NAME instead of a per-sample Control column for one shared control.
# Outputs: <outprefix>_gene_JACKS_results.txt, _gene_std_, _grna_ (and _gene_pval_ with --ctrl-genes
#          AND --n-pseudo > 0), _JACKS_results_full.pickle.
import argparse
import random
from jacks.jacks_io import runJACKS

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('counts')
ap.add_argument('replicatemap')
ap.add_argument('guidemap')
ap.add_argument('--outprefix', default='jacks_out')
ap.add_argument('--rep-hdr', default='Replicate')
ap.add_argument('--sample-hdr', default='Sample')
ap.add_argument('--ctrl-sample-hdr', default='Control', help='per-sample control column')
ap.add_argument('--common-ctrl-sample', default=None, help='one shared control sample name instead')
ap.add_argument('--sgrna-hdr', default='sgRNA')
ap.add_argument('--gene-hdr', default='Gene')
ap.add_argument('--ctrl-genes', default=None, help='negative-control gene list; required for p-values')
ap.add_argument('--n-pseudo', type=int, default=0, help='pseudo-genes for p-values; the Python default is 0 (CLI: 2000)')
ap.add_argument('--seed', type=int, default=None, help='seeds random; for byte-identical p-values also set PYTHONHASHSEED in the environment')
ap.add_argument('--apply-w-hp', action='store_true',
                help='hierarchical gene-effect prior; the tool advises caution, off by default')
a = ap.parse_args()

if a.seed is not None:
    random.seed(a.seed)                       # only the pseudo-gene p-values are random

kwargs = dict(
    countfile=a.counts,
    replicatefile=a.replicatemap,
    guidemappingfile=a.guidemap,
    rep_hdr=a.rep_hdr,
    sample_hdr=a.sample_hdr,
    sgrna_hdr=a.sgrna_hdr,
    gene_hdr=a.gene_hdr,
    outprefix=a.outprefix,
    apply_w_hp=a.apply_w_hp,                  # False is the default and the tool's recommendation
)
if a.common_ctrl_sample:
    kwargs['common_ctrl_sample'] = a.common_ctrl_sample
else:
    kwargs['ctrl_sample_hdr'] = a.ctrl_sample_hdr
if a.ctrl_genes:
    kwargs['ctrl_genes'] = a.ctrl_genes
    kwargs['n_pseudo'] = a.n_pseudo
runJACKS(**kwargs)
