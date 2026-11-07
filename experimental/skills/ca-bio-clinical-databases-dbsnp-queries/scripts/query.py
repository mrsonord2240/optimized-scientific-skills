"""Research-only public dbSNP example queries."""
import argparse
import json
from pathlib import Path
from dbsnp_lookup import resolve_merge_chain, vcf_to_canonical_spdi, alfa_population_frequencies

p = argparse.ArgumentParser(description=__doc__)
sub = p.add_subparsers(dest="mode", required=True)
r = sub.add_parser("resolve")
r.add_argument("rsid")
r.add_argument("output")
v = sub.add_parser("vcf")
v.add_argument("chrom")
v.add_argument("pos", type=int)
v.add_argument("ref")
v.add_argument("alt")
v.add_argument("output")
a = sub.add_parser("alfa")
a.add_argument("rsid")
a.add_argument("output")
args = p.parse_args()
if args.mode == "resolve":
    result = resolve_merge_chain(args.rsid)
elif args.mode == "vcf":
    result = vcf_to_canonical_spdi(args.chrom, args.pos, args.ref, args.alt)
else:
    result = alfa_population_frequencies(args.rsid)
Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
