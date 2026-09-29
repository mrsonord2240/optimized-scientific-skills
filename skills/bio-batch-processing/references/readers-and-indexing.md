# Reader selection and indexing

Read this reference when choosing a parser, planning repeated access, handling
compressed inputs, or diagnosing reader-specific failures.

## Dependencies

Biopython is required for the core recipes. pysam and pyfastx are optional
alternatives for large FASTQ or indexed compressed FASTX:

```bash
pip install biopython
pip install pysam pyfastx
```

The source examples target Biopython 1.83+, pysam 0.22+, and pyfastx 2.0+.

## Why streaming is the default

`SeqIO.parse()` returns a generator that holds one `SeqRecord` at a time.
`list(SeqIO.parse(...))` materializes every record and can exhaust memory when
applied to large files or a directory. Use a list only when the input is known
to be small and multiple in-memory passes are required.

For random access, `SeqIO.index_db()` builds a persistent SQLite index over one
or more files. It reparses only the selected record on lookup and avoids loading
the corpus into memory.

At tens of millions of reads, full `SeqRecord` construction can dominate runtime.
When the task is plain iteration, choose a thinner reader.

## Reader matrix

| Reader | Per-record object | Random access | Best fit |
|---|---|---|---|
| `Bio.SeqIO.parse` | Full `SeqRecord` | No; one-pass generator | Small or medium data that needs the Biopython record API |
| `Bio.SeqIO.index_db` | Reparsed `SeqRecord` on access | Yes; persistent, on-disk, multi-file | Repeated lookup across many or large files |
| `pysam.FastxFile` | Thin entry with name, sequence, comment, and quality | No; linear access, including sequential gzip | Fast linear iteration over huge FASTQ |
| `pyfastx` | Tuple or object backed by a SQLite index | Yes; plain or gzipped FASTA/FASTQ | Indexed reuse of compressed or uncompressed FASTX |

`pysam.FastxFile.get_quality_array()` always subtracts 33. It is correct for
Phred+33 data, not legacy Phred+64 or Solexa encodings. For those encodings use
`SeqIO` with the explicit format variant.

pyfastx creates persistent `.fxi` or `.fqi` SQLite indexes and can randomly
access plain or gzipped FASTA/FASTQ without re-compressing them as BGZF.

### Runnable pyfastx gzip-index example

From the Skill directory, build or reuse the persistent `.fxi`, retrieve a
checked id, and print JSON describing the record count and selected length:

```bash
python scripts/pyfastx_index.py reference.fasta.gz --record-id seq_00042
```

The script fails if the FASTA is empty, the requested id is absent, or the
expected sidecar was not created. Current pyfastx releases do not expose a
public close method on every indexed object; the script calls it when present
and otherwise releases the object at process/function exit. The `.fxi` is a
deliberate reusable artifact—delete it only when invalidating that index.

## Persistent random access across files

Use one index for the file set:

```python
from pathlib import Path
from Bio import SeqIO

def stable_key(path):
    relative = path.relative_to("data").as_posix()
    return relative.casefold(), relative

files = [str(path) for path in sorted(Path("data").glob("*.fasta"), key=stable_key)]
records = SeqIO.index_db("combined.idx", files, "fasta")

print(len(records))
record = records["seq_00042"]
records.close()
```

The index persists. A later process can reopen it without passing the file list:

```python
from Bio import SeqIO

records = SeqIO.index_db("combined.idx")
try:
    print(records["seq_00042"])
finally:
    records.close()
```

Record ids must be unique across every indexed file. A collision raises
`ValueError: Duplicate key`. Prefix ids by source file or provide an appropriate
`key_function` when the ids are not globally unique.

`SeqIO.index_db()` supports uncompressed and BGZF-compressed inputs. Plain gzip
is not seekable for this index. Recompress with `bgzip`, or use pyfastx when its
gzip indexing semantics match the task.

## Huge FASTQ linear count

When only thin FASTQ fields are needed:

```python
import pysam

with pysam.FastxFile("reads.fastq.gz") as handle:
    count = sum(1 for _ in handle)
```

## Diagnostic table

| Symptom | Likely cause | Response |
|---|---|---|
| `MemoryError` or killed process | Parsed records were materialized | Stream the generator and eliminate whole-file lists |
| Linear job is slow on tens of millions of reads | Full `SeqRecord` construction is unnecessary overhead | Use `pysam.FastxFile`, or pyfastx for indexed access |
| `ValueError: Duplicate key` | The same id occurs in more than one indexed input | Make ids globally unique or use a suitable key function |
| A second pass yields no records | The generator was exhausted | Create a new parser or use an index |
| `index_db` fails on `.gz` | Input is plain gzip rather than BGZF | Recompress as BGZF or choose pyfastx |
| Converted output lacks annotations | The target format cannot represent them | Retain an annotated format or extract qualifiers before conversion |
