"""Public-example ClinVar or local snapshot lookup; not clinical interpretation."""
import argparse
import json
from pathlib import Path
from clinvar_query import clinvar_summary, car_record, lookup_local_vcf

p = argparse.ArgumentParser(description=__doc__)
sub = p.add_subparsers(dest="mode", required=True)
s = sub.add_parser("summary")
s.add_argument("variation_id")
s.add_argument("output")
a = sub.add_parser("allele")
a.add_argument("hgvs")
a.add_argument("output")
v = sub.add_parser("local")
v.add_argument("vcf")
v.add_argument("chrom")
v.add_argument("pos", type=int)
v.add_argument("ref")
v.add_argument("alt")
v.add_argument("output")
args = p.parse_args()
if args.mode == "summary":
    result = clinvar_summary(args.variation_id)
elif args.mode == "allele":
    result = car_record(args.hgvs)
else:
    result = lookup_local_vcf(args.vcf, args.chrom, args.pos, args.ref, args.alt)
Path(args.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
