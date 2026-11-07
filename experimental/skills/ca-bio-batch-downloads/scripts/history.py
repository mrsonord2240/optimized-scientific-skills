#!/usr/bin/env python3
"""Run source checkpointed Entrez retrieval for user-specified inputs."""
import argparse
import os
from pathlib import Path
from Bio import Entrez
from robust_download import checkpointed_download

def main():
    p = argparse.ArgumentParser()
    p.add_argument('db')
    p.add_argument('query')
    p.add_argument('output')
    p.add_argument('checkpoint')
    p.add_argument('--email', required=True)
    p.add_argument('--rettype', default='fasta')
    p.add_argument('--batch-size', type=int, default=500)
    args = p.parse_args()
    if args.batch_size < 1:
        p.error('batch size must be positive')
    if Path(args.output).resolve() == Path(args.checkpoint).resolve():
        p.error('output and checkpoint must differ')
    if Path(args.output).exists() and not Path(args.checkpoint).exists():
        p.error('output exists without checkpoint; choose a fresh path')
    if Path(args.checkpoint).exists() and not Path(args.output).exists():
        p.error('checkpoint exists without its output')
    Entrez.email = args.email
    Entrez.api_key = os.environ.get('NCBI_API_KEY')
    checkpointed_download(args.db, args.query, args.output, args.checkpoint,
                          rettype=args.rettype, batch_size=args.batch_size)

if __name__ == '__main__':
    main()
