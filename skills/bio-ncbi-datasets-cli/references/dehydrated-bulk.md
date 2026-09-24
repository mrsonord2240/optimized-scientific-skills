# Dehydrated bulk downloads and checksum verification

Read when pulling many genomes (hundreds or more), transferring with aria2c, or checking downloaded files.

## When to use --dehydrated for cloud workflows

The "dehydrated" mode separates data discovery from data transfer:

1. **Discover**: `datasets download genome taxon human --reference --dehydrated --filename human.zip` (fast; ~MB).
2. **Inspect**: `unzip -p human.zip ncbi_dataset/fetch.txt` -- a TSV of all URLs to pull.
3. **Pull**: either `datasets rehydrate --directory ./human/`, or `aria2c` for parallel pull.
   `fetch.txt` is 3 tab-separated columns (`<url>`, a `0` placeholder, `<path relative to
   ncbi_dataset/>`), not aria2c's input format -- convert it first (see the bulk pattern below).
   After an aria2c pull, size-check the files (see Checksum verification): `rehydrate` will not.

This is essential for HPC / cloud pipelines where inspection of the pending transfer is needed before committing the I/O.

## Checksum verification

`datasets download` validates the downloaded zip's checksum by default (`--fast-zip-validation`
skips it). This replaces the `md5sum -c` step that assembly_summary.txt-based scraping needed.

**`datasets rehydrate` does not verify anything already on disk** (checked on 18.37.0, 2026-09-21).
It downloads only files missing from `ncbi_dataset/data/`; a file that exists at the expected path
counts as "already rehydrated" whatever its content. After an `aria2c` pull that is the dangerous
case: a throttled or blocked transfer can write an HTML error page (a few KB) at the correct path,
`aria2c` reports success, and rehydrate says `All N files already rehydrated`. Reproduced by
overwriting a rehydrated `.fna` with 4 bytes -- rehydrate left it untouched.

Check sizes yourself against the byte lengths the CLI recorded in `dataset_catalog.json`
(`uncompressedLengthBytes`), delete any
mismatching file, and rehydrate again to re-fetch just those:

```bash
python3 - ncbi_dataset/data <<'PY'    # arg: the dehydrated package's ncbi_dataset/data dir
import json, os, sys
root = sys.argv[1]
for asm in json.load(open(os.path.join(root, 'dataset_catalog.json')))['assemblies']:
    for f in asm['files']:
        p = root + '/' + f['filePath']
        if not os.path.exists(p) or os.path.getsize(p) != int(f['uncompressedLengthBytes']):
            print(p)                  # bad or missing: rm it, then `datasets rehydrate --directory <pkg>`
PY
```

`examples/bulk_dehydrated.sh` runs this check and the delete-and-rehydrate retry automatically.

### Bulk download all reference bacterial genomes

**Goal:** Pull every RefSeq reference bacterial assembly with annotation.

**Approach:** `--dehydrated` first for inspection; rehydrate with parallel pull.

`examples/bulk_dehydrated.sh` runs all three steps (checked on NCBI Datasets CLI 18.37.0, 2026-09-21):
dehydrated discovery, the `fetch.txt` -> aria2c input conversion (path is column 3, `--dir` is
`ncbi_dataset/`) with an `out=0` sanity check, then the size check and delete-and-rehydrate retry.

```bash
examples/bulk_dehydrated.sh Bacteria bact_refs.zip bact_refs 8    # taxon, zip, dest dir, threads
```
