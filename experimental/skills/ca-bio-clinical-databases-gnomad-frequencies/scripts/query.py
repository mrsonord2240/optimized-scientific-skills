"""Public example gnomAD lookup; not for clinical interpretation."""
import argparse
import json
from pathlib import Path
from gnomad_query import query_variant_v4, grpmax_faf95, query_gene_constraint_v4

p = argparse.ArgumentParser(description=__doc__)
sub = p.add_subparsers(dest="mode", required=True)
v = sub.add_parser("variant")
v.add_argument("chrom")
v.add_argument("pos", type=int)
v.add_argument("ref")
v.add_argument("alt")
v.add_argument("output")
g = sub.add_parser("gene")
g.add_argument("symbol")
g.add_argument("output")
args = p.parse_args()
if args.mode == "variant":
    raw = query_variant_v4(args.chrom, args.pos, args.ref, args.alt)
    result = {"variant": raw, "faf95": grpmax_faf95(raw)}
else:
    result = query_gene_constraint_v4(args.symbol)
Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
