#!/usr/bin/env python3
"""Build/reuse a pyfastx FASTA index and perform checked random access."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyfastx


def inspect_fasta(input_path, record_id=None):
    input_path = Path(input_path)
    fasta = pyfastx.Fasta(str(input_path), build_index=True, full_name=False)
    try:
        keys = fasta.keys()
        if len(keys) == 0:
            raise ValueError(f"FASTA contains no records: {input_path}")
        selected_id = record_id or keys[0]
        record = fasta[selected_id]
        result = {
            "input": str(input_path),
            "index": f"{input_path}.fxi",
            "records": len(fasta),
            "selected_id": selected_id,
            "selected_length": len(record),
        }
    finally:
        close = getattr(fasta, "close", None)
        if callable(close):
            close()

    if not Path(result["index"]).is_file():
        raise RuntimeError(f"pyfastx did not create the expected index: {result['index']}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="plain or gzip-compressed FASTA")
    parser.add_argument("--record-id", help="record to retrieve; default: first indexed id")
    args = parser.parse_args()
    print(json.dumps(inspect_fasta(args.input, args.record_id), sort_keys=True))


if __name__ == "__main__":
    main()
