"""Research-only public myvariant.info batch and gene queries."""
import argparse
import json
from pathlib import Path
from query_myvariant import annotate_variants, find_pathogenic_in_gene

p = argparse.ArgumentParser(description=__doc__)
sub = p.add_subparsers(dest="mode", required=True)
b = sub.add_parser("batch")
b.add_argument("input")
b.add_argument("output")
b.add_argument("versions")
g = sub.add_parser("gene")
g.add_argument("symbol")
g.add_argument("output")
g.add_argument("limit", type=int, nargs="?", default=500)
g.add_argument("min_star", type=int, nargs="?", default=2)
args = p.parse_args()
if args.mode == "batch":
    ids = [x.strip() for x in Path(args.input).read_text(encoding="utf-8").splitlines() if x.strip()]
    if not ids:
        p.error("input must contain at least one identifier")
    frame, versions = annotate_variants(ids)
    Path(args.output).write_text(frame.to_json(orient="records", indent=2) + "\n", encoding="utf-8")
    Path(args.versions).write_text(json.dumps(versions, indent=2) + "\n", encoding="utf-8")
else:
    result = find_pathogenic_in_gene(args.symbol, args.limit, args.min_star)
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
