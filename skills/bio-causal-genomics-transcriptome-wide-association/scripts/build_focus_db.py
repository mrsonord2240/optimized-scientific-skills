#!/usr/bin/env python
"""Build a FOCUS (pyfocus) weight database directly from a weight table.

Bypasses `focus import ... fusion`, which needs mygene + rpy2 (rpy2 needs R built as a shared
library) and otherwise leaves an empty DB with only an ERROR log line. `focus finemap` reads the
result as a normal DB.

Input  : panel.tsv, tab-separated, one row per gene-SNP weight, columns
         gene chrom txstart txstop snp pos a1 a0 weight
         (optional columns cv_r2, cv_r2_pval carry each gene's real CV R2 and p; default 0.1 and 1e-4)
Output : <out.db> (pyfocus sqlite schema, one RefPanel + one model per gene)
Usage  : python build_focus_db.py panel.tsv custom_focus.db [--ref-name custom_panel] [--tissue Whole_Blood]
Tested : pyfocus 0.802, pandas 2.1.4, SQLAlchemy 2.0.54 (see SKILL.md Tool Install Notes for the pins/patch).
"""
import argparse

import pandas as pd
from pyfocus.models.db import load_db, RefPanel, build_model

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("panel_tsv")
ap.add_argument("out_db")
ap.add_argument("--ref-name", default="custom_panel")
ap.add_argument("--tissue", default="Whole_Blood")
ap.add_argument("--assay", default="rnaseq")
args = ap.parse_args()

panel = pd.read_csv(args.panel_tsv, sep="\t")
ssn = load_db(args.out_db)
ref = RefPanel(ref_name=args.ref_name, tissue=args.tissue, assay=args.assay)
for gene, g in panel.groupby("gene"):
    snps = g.reset_index(drop=True)
    snps["chrom"] = snps["chrom"].astype(str)
    r2 = float(g["cv_r2"].iloc[0]) if "cv_r2" in g else 0.1        # use each gene's real CV R2 and p
    pv = float(g["cv_r2_pval"].iloc[0]) if "cv_r2_pval" in g else 1e-4
    ssn.add(build_model(
        gene_info=dict(geneid=gene, txid=gene, name=gene, type="protein_coding",
                       chrom=str(g.chrom.iloc[0]), txstart=int(g.txstart.iloc[0]), txstop=int(g.txstop.iloc[0])),
        snp_info=snps, db_ref_panel=ref, weights=g.weight.tolist(), ses=None,
        attrs={"cv.R2": r2, "cv.R2.pval": pv}, method="top1"))
ssn.add(ref)
ssn.commit()
print(f"Wrote {args.out_db}: {panel.gene.nunique()} genes, {len(panel)} weights")
