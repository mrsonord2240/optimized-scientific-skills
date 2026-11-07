#!/usr/bin/env python3
"""File-oriented entry points for the source cfDNA preprocessing functions."""
import argparse
import json
from pathlib import Path
from preprocess_cfdna import insert_size_qc, preprocess_cfdna


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    consensus = sub.add_parser("consensus")
    consensus.add_argument("input_bam")
    consensus.add_argument("output_bam")
    consensus.add_argument("reference")
    consensus.add_argument("--read-structures", nargs=2, default=("6M11S+T", "6M11S+T"))
    consensus.add_argument("--tags", nargs="+", default=("ZA", "ZB"))
    consensus.add_argument("--threads", type=int, default=8)
    consensus.add_argument("--duplex", action="store_true")
    qc = sub.add_parser("qc")
    qc.add_argument("input_bam")
    qc.add_argument("output_json")
    qc.add_argument("--max-size", type=int, default=600)
    args = parser.parse_args()
    if args.command == "consensus":
        preprocess_cfdna(args.input_bam, args.output_bam, args.reference,
                         read_structure=args.read_structures, duplex=args.duplex,
                         threads=args.threads, molecular_index_tags=args.tags)
    else:
        if Path(args.input_bam).resolve() == Path(args.output_json).resolve():
            parser.error("output JSON must differ from input BAM")
        result = insert_size_qc(args.input_bam, max_size=args.max_size)
        Path(args.output_json).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
