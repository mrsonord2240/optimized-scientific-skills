# Batch operation recipes

Read only the section for the requested operation. Run commands from the Skill
directory so `scripts.batch_process` is importable. The shared functions are
covered by `python -B tests/test_batch_process.py` and require Biopython 1.83+.

## Deterministic discovery

Every multi-file operation must use the same discovery order. `stable_paths`
sorts relative POSIX paths by `(casefolded path, exact path)`, so filesystem
creation order cannot change merge bytes, summary rows, conversions, indexes,
or parallel result order.

```python
from scripts.batch_process import stable_paths

files = stable_paths("data", "*.fasta")
recursive_genbank = stable_paths("data", "*.gb", recursive=True)
```

## Count and merge

### Count records across files

Count by consuming each parser rather than building a record list:

```python
from scripts.batch_process import count_streaming, stable_paths

for fasta_file in stable_paths("data", "*.fasta"):
    print(f"{fasta_file.name}: {count_streaming(fasta_file, 'fasta')} sequences")
```

### Merge without materializing records

```python
from Bio import SeqIO
from scripts.batch_process import stable_paths

def all_records(directory, pattern, format):
    for filepath in stable_paths(directory, pattern):
        with filepath.open("r", encoding="utf-8") as handle:
            yield from SeqIO.parse(handle, format)

count = SeqIO.write(
    all_records("data", "*.fasta", "fasta"),
    "merged.fasta",
    "fasta",
)
print(f"Merged {count} records")
```

### Merge with source tracking

```python
from Bio import SeqIO
from scripts.batch_process import stable_paths

def records_with_source(directory, pattern, format):
    for filepath in stable_paths(directory, pattern):
        with filepath.open("r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, format):
                record.description = f"{record.description} [source={filepath.name}]"
                yield record

SeqIO.write(
    records_with_source("data", "*.fasta", "fasta"),
    "merged_tracked.fasta",
    "fasta",
)
```

A streamed FASTA merge permits duplicate ids. A later `index_db()` or
`to_dict()` over that output does not; decide whether to prefix ids before
creating a merged artifact intended for indexed access.

## Split files

### Split by record count

`split_by_count` rejects Boolean, non-integer, zero, and negative chunk sizes
before it opens the input or creates output. A positive chunk retains at most
`records_per_file` records.

```python
from scripts.batch_process import split_by_count

outputs = split_by_count("large.fasta", "fasta", 1000, "split")
print(f"Wrote {len(outputs)} chunks")
```

### Split by sequence-id prefix

Treat record ids as untrusted data. This function accepts only portable ASCII
prefixes, rejects path separators and reserved filenames, verifies target
containment, and creates each output exclusively. It keeps at most 64 handles
open by default and writes `split-prefix-manifest.json` with `complete` or
`partial` status, per-file counts, and any error. Existing outputs or manifests
are never overwritten; choose a new output directory for a new attempt.

```python
from scripts.batch_process import split_by_prefix

result = split_by_prefix(
    "input.fasta",
    "fasta",
    "split-by-prefix",
    max_open_files=64,
)
assert result["status"] == "complete"
print(result["records_written"], len(result["files"]))
```

The prefix is the part of `record.id` before the first underscore. If valid ids
need a different grouping rule, transform them explicitly to a reviewed,
collision-free portable token before invoking the split.

## Batch conversion

`SeqIO.convert()` streams internally. Discover inputs deterministically and
create the destination directory before processing:

```python
from pathlib import Path
from Bio import SeqIO
from scripts.batch_process import stable_paths

destination = Path("fasta")
destination.mkdir(parents=True, exist_ok=True)
for gb_file in stable_paths("genbank", "*.gb"):
    fasta_file = destination / gb_file.with_suffix(".fasta").name
    count = SeqIO.convert(str(gb_file), "genbank", str(fasta_file), "fasta")
    print(f"{gb_file.name} -> {fasta_file.name}: {count} records")
```

GenBank-to-FASTA conversion drops features, annotations, and qualifiers because
FASTA stores only id, description, and sequence. Keep an annotated format or
extract the required metadata first.

## Parallel processing

For CPU-bound per-file work, distribute whole files across processes. The
worker must remain importable at module scope, while Pool construction and the
program call stay behind the main-module guard required by spawn platforms.

```python
from scripts.batch_process import process_files_parallel, stable_paths

def main():
    files = stable_paths("data", "*.fasta")
    results = process_files_parallel(files, workers=4)
    print(results)

if __name__ == "__main__":
    main()
```

The file-path list is materialized, not the sequence records. `results` retains
the stable input order under both fork and spawn. Use
`concurrent.futures.ThreadPoolExecutor` only for genuinely I/O-bound work; the
Python GIL limits threads for CPU-bound parsing.

## Summary statistics

`summarize_files` streams each file once, uses the stable input order, and
writes the explicit columns `file,sequences,total_bp,min_len,max_len,avg_len`.
No matches produce a valid header-only CSV and an empty result list.

```python
from scripts.batch_process import stable_paths, summarize_files

files = stable_paths("data", "*.fasta")
summaries = summarize_files(files, "fasta", "summary.csv")
print(f"Summarized {len(summaries)} files")
```
