#!/usr/bin/env python3
"""File-based entry points for streaming batch operations."""
import argparse
import json
from pathlib import Path
from Bio import SeqIO
from batch_process import stable_paths, count_streaming, split_by_count, split_by_prefix, summarize_files, process_files_parallel

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('operation', choices=['count', 'merge', 'summary', 'convert', 'split', 'prefix', 'parallel'])
    p.add_argument('input')
    p.add_argument('pattern', nargs='?', default='*.fasta')
    p.add_argument('--format', default='fasta')
    p.add_argument('--output')
    p.add_argument('--output-format')
    p.add_argument('--chunk-size', type=int)
    p.add_argument('--workers', type=int, default=4)
    a = p.parse_args()
    if a.operation in ('merge', 'summary', 'convert', 'split', 'prefix') and not a.output:
        p.error('--output is required for this operation')
    if a.operation == 'split':
        if a.chunk_size is None:
            p.error('--chunk-size is required for split')
        print(json.dumps(split_by_count(a.input, a.format, a.chunk_size, a.output)))
        return
    if a.operation == 'prefix':
        print(json.dumps(split_by_prefix(a.input, a.format, a.output)))
        return
    files = stable_paths(a.input, a.pattern)
    if a.operation == 'count':
        print(json.dumps([{'file': str(f), 'count': count_streaming(f, a.format)} for f in files]))
    elif a.operation == 'parallel':
        print(json.dumps(process_files_parallel(files, workers=a.workers)))
    elif a.operation == 'summary':
        summarize_files(files, a.format, a.output)
    elif a.operation == 'merge':
        def records():
            for f in files:
                with f.open('r', encoding='utf-8') as h:
                    yield from SeqIO.parse(h, a.format)
        print(SeqIO.write(records(), a.output, a.format))
    elif a.operation == 'convert':
        if not a.output_format:
            p.error('--output-format is required for convert')
        dest = Path(a.output)
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            target = dest / (f.stem + '.' + a.output_format)
            print(f'{f} -> {target}: {SeqIO.convert(str(f), a.format, str(target), a.output_format)}')

if __name__ == '__main__':
    main()
