#!/usr/bin/env python3
"""Write first-HSP summaries from a single saved BLAST XML result."""
import csv
import sys
from Bio.Blast import NCBIXML

if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('usage: parse_results.py input.xml output.tsv')
    with open(sys.argv[1], encoding='utf-8') as handle:
        record = NCBIXML.read(handle)
    with open(sys.argv[2], 'w', encoding='utf-8', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['accession', 'title', 'evalue', 'bits', 'identity', 'source_coverage', 'query_start', 'query_end'])
        for alignment in record.alignments:
            if not alignment.hsps:
                continue
            h = alignment.hsps[0]
            writer.writerow([alignment.accession, alignment.title, h.expect, h.bits,
                             h.identities / h.align_length if h.align_length else '',
                             h.align_length / record.query_length if record.query_length else '',
                             h.query_start, h.query_end])
