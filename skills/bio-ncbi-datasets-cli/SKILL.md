---
name: bio-ncbi-datasets-cli
description: Download genome assemblies, gene records, and ortholog data from NCBI using the modern Datasets v2 CLI (replaces assembly_summary.txt scraping and many EFetch workflows). Use when bulk-pulling genome assemblies, gene metadata across species, ortholog sets, or BLAST databases; when E-utilities are too slow for genome-scale work; or when download-time zip checksum validation, parallel download, and clean accession-driven retrieval are required. Encodes the JSON-lines output format, dataformat conversion, --dehydrated for cloud workflows, and when Datasets is/isn't the right tool.
tool_type: cli
primary_tool: NCBI Datasets CLI
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples checked on NCBI Datasets CLI 18.37.0 / dataformat 18.37.0 (2026-09-19). Earlier
docs pinned this Skill to 16.0+ (2024); the 16 -> 18 jump renamed the `dataformat --fields` catalog
and changed `--ortholog` from a boolean flag to `--ortholog strings` (see Code patterns below).

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `datasets --version` (reliable). `dataformat --version` is broken on 18.37.0 -- it prints the
  literal string `undefined` regardless of the real version. Confirm the `dataformat` build instead
  via `dataformat --help`'s banner, or by checking it shipped alongside a known-good `datasets`
  binary (they are always released as a matched pair).
- Subcommand help: `datasets <subcommand> --help`

If a subcommand or flag is unrecognized, run `datasets --help` and adapt. The CLI is under active
development; major releases have added subcommands and renamed flags and fields.

# NCBI Datasets CLI

**"Pull genome / gene / ortholog data from NCBI in 2026"** -> The Datasets v2 CLI (launched 2023) is the official, supported bulk endpoint for genome and gene-centric data. It replaces the prior best-practice of scraping `assembly_summary.txt` + parallel FTP + manual checksum verification. For genome-scale data, it is strictly better than E-utilities (EFetch).

The CLI is not the right answer for everything. PubMed, SRA reads, and custom Entrez queries still belong to E-utilities. The defection rule: **if the question is about genome assemblies, gene records, or pre-computed orthologs, use Datasets; otherwise stay with E-utilities**.

- CLI: `datasets download genome accession GCF_...`
- CLI: `datasets summary gene symbol BRCA1 --taxon human`
- Python: `subprocess` wrapper; Python client `ncbi-datasets-pylib` (experimental as of 2024)

## Installation

```bash
# conda
conda install -c conda-forge ncbi-datasets-cli

# Or direct download (Linux, macOS, Windows binaries)
curl -O https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets

datasets --version    # this doc was checked against 18.37.0
dataformat --version  # bundled companion tool; broken on 18.37.0, see Version Compatibility
```

## What's in scope (use Datasets) vs out of scope (use E-utilities or other tools)

| Question | Datasets | Use instead |
|---|---|---|
| Genome assembly download | yes | — |
| All reference genomes for a taxon | yes | — |
| Gene record metadata (multi-species) | yes (see `references/gene-orthologs.md`) | — |
| Ortholog data for a gene | yes (`datasets summary gene ... --ortholog <taxon|all>`; see `references/gene-orthologs.md`) | OrthoDB / Compara for tree-aware orthology |
| Virus data (assemblies, metadata) | yes (`datasets download virus`; see `references/virus-genomes.md`) | — |
| Annotation files (GFF3, GTF) for a genome | yes | — |
| Protein records (curated, with cross-refs) | partial | UniProt REST for richer annotation |
| PubMed | no | `entrez-search` / `entrez-fetch` |
| SRA reads | no | `sra-data` |
| BLAST | no | `blast-searches` / `local-blast` |
| Custom Entrez queries | no | `entrez-search` |
| Pre-computed alignments (Compara) | no | `ensembl-rest` |

## Subcommand taxonomy

| Subcommand | Purpose | Example |
|---|---|---|
| `datasets summary genome` | Metadata only; JSON output | `datasets summary genome accession GCF_000001405.40` |
| `datasets download genome` | Download data files | `datasets download genome accession GCF_...` |
| `datasets summary gene` | Gene record metadata | `datasets summary gene symbol BRCA1 --taxon human` |
| `datasets download gene` | Download gene products | `datasets download gene symbol BRCA1 --taxon human` |
| `datasets summary taxonomy` | Taxonomy info | `datasets summary taxonomy taxon human` |
| `datasets download virus` | Virus assemblies/proteins | `datasets download virus genome taxon SARS-CoV-2` |
| `dataformat tsv` / `dataformat excel` | Convert JSON-lines to tabular | `dataformat tsv gene-summary` |

`datasets summary` always returns JSON-lines on stdout (one object per record). `datasets download` produces a `.zip` (default) or a "dehydrated" stub for cloud workflows.

## Key parameters (download)

| Flag | Effect |
|---|---|
| `--filename out.zip` | Where to write the archive |
| `--include genome,gff3,gtf,protein,cds,rna,seq-report` | Which file types to include |
| `--reference` | Restrict to reference assemblies only (one per species) |
| `--annotated` | Restrict to annotated assemblies |
| `--assembly-source RefSeq` / `GenBank` / `all` | Database source |
| `--assembly-level chromosome,complete` | Assembly quality level |
| `--released-after 2024-01-01` | Date filter |
| `--dehydrated` | Skip data; download just stubs + URL list (for parallel pull) |
| `--api-key XXX` | Optional API key (raises rate limit) |
| `--no-progressbar` | For non-interactive use |

For very large pulls (1000+ genomes), `--dehydrated` is the right choice: download the metadata stubs first, then run `datasets rehydrate` later or pull URLs in parallel from the manifest. The workflow, the aria2c conversion of `fetch.txt` and the post-transfer size check (`rehydrate` does not verify existing files) are in `references/dehydrated-bulk.md`. `download` itself validates the zip checksum (`--fast-zip-validation` skips it).

## JSON-lines output + dataformat

`datasets summary` returns JSON-lines (one JSON object per line) on stdout. Pipe through `dataformat tsv` for tabular:

```bash
datasets summary genome taxon "Escherichia coli" --reference --as-json-lines \
  | dataformat tsv genome --fields accession,organism-name,assminfo-level,assmstats-scaffold-n50 \
  > ecoli_refs.tsv
```

`dataformat` subcommands match summary types: `genome`, `gene`, `virus-genome`, etc. The `--fields`
list is documented per type via `dataformat tsv <type> --help` -- **always check it live**, do not
reuse a list from an older doc. On 18.37.0 the genome catalog uses `assminfo-*` / `assmstats-*`
prefixes (e.g. `assminfo-level`, `assmstats-scaffold-n50`, `assmstats-contig-n50`,
`assmstats-total-sequence-len`) and the gene catalog uses `tax-name` (not `taxname`); there is no
`nomenclature-authority-symbol` field on this build.

## Reference Files

| File | Read when |
|---|---|
| `references/dehydrated-bulk.md` | Pulling hundreds of genomes or more, transferring with aria2c or `datasets rehydrate`, or verifying files after a transfer |
| `references/gene-orthologs.md` | A gene request spans more than one species, or asks for orthologs (`--ortholog`) |
| `references/virus-genomes.md` | Virus assemblies, metadata or proteins (`datasets download virus`) |

## Code patterns

### Download a single reference genome

**Goal:** Get human reference assembly with genome + GTF + protein + CDS.

**Approach:** `datasets download genome accession ... --include ...`.

**Reference (NCBI Datasets CLI 18.37.0, checked 2026-09-19):**
```bash
#!/bin/bash
# Reference: NCBI Datasets CLI 18.37.0 (checked 2026-09-19) | Verify API if version differs

datasets download genome accession GCF_000001405.40 \
    --include genome,gff3,gtf,protein,cds,seq-report \
    --filename human_grch38.zip

unzip -q human_grch38.zip -d human_grch38/
ls -lh human_grch38/ncbi_dataset/data/GCF_000001405.40/
```

### Filter assemblies by quality and date

```bash
datasets summary genome taxon "Salmonella enterica" \
    --assembly-level chromosome,complete \
    --released-after 2024-01-01 \
    --as-json-lines \
  | dataformat tsv genome --fields accession,organism-name,assminfo-level,assmstats-scaffold-n50,assminfo-release-date \
  > sal_2024.tsv
```

### Python wrapper

`scripts/datasets_wrapper.py` (checked on NCBI Datasets CLI 18.37.0, 2026-09-19) wraps `datasets summary`
(JSON-lines parsed into dicts; keys are snake_case, e.g. `assembly_stats.contig_n50`) and
`datasets download` (returns the zip path). Import `datasets_summary` / `datasets_download`, or run:

```bash
python scripts/datasets_wrapper.py --taxon "Escherichia coli" --accession GCF_000005845.2 --out ecoli_k12.zip --include genome,gff3,protein
```

### Comparison vs E-utilities

```python
# E-utilities path: ESearch in assembly db -> ESummary -> manual FTP pull
#   ~30 API calls + manual md5 + serial download
# Datasets path:
#   datasets download genome accession GCF_...  # one command, zip checksum validated
```

For genome workflows, Datasets is 5-50x faster than the equivalent E-utilities pipeline and far more reliable.

## Failure modes

### Choosing Datasets for the wrong question
- **Trigger:** Trying to pull raw SRA reads via Datasets.
- **Mechanism:** Datasets covers genome/gene/ortholog, not raw reads.
- **Symptom:** Subcommand not found or empty result.
- **Fix:** Use `sra-data` skill (prefetch/fasterq-dump) for raw reads.

### `--reference` filter loses too much
- **Trigger:** Bulk pull of "all assemblies for a species"; `--reference` returns one per species.
- **Mechanism:** Reference subset is the canonical single representative.
- **Symptom:** Far fewer assemblies than expected for a species with hundreds of submissions.
- **Fix:** Drop `--reference` for full set; add `--assembly-level chromosome,complete` for quality filter instead.

### Dehydrated workflow forgotten
- **Trigger:** 1000-genome pull without `--dehydrated`.
- **Mechanism:** Datasets downloads serially within one ZIP; can take hours.
- **Symptom:** Slow; no parallelism; one giant ZIP.
- **Fix:** Use `--dehydrated` + aria2c with `--max-concurrent-downloads`.

### dataformat field name guessing
- **Trigger:** `dataformat tsv genome --fields foo,bar` with invented field names.
- **Mechanism:** Field names are constrained per summary type.
- **Symptom:** "Unknown field" error.
- **Fix:** `dataformat tsv genome --help` lists valid field names; pull JSON-lines and inspect with `jq` to discover fields.

### Old assembly_summary.txt-based scripts still in use
- **Trigger:** Legacy pipeline scraping `https://ftp.ncbi.nlm.nih.gov/genomes/all/refseq/...`.
- **Mechanism:** Pre-2023 best practice; FTP listing parsing is fragile.
- **Symptom:** Slow; brittle; no checksums; broken when NCBI restructures FTP.
- **Fix:** Switch to Datasets CLI; the FTP path still works but Datasets is the supported modern path.

### API key not used for high-volume
- **Trigger:** 1000+ summary calls in a loop without `--api-key`.
- **Mechanism:** NCBI rate-limits unauthenticated bulk traffic.
- **Symptom:** Throttling; slow downloads.
- **Fix:** Pass `--api-key YOUR_KEY` to bulk commands; obtain from `https://www.ncbi.nlm.nih.gov/account/settings/`.

### CLI version drift
- **Trigger:** Running docs written for one CLI version against a different installed version (this
  doc was last checked against 18.37.0; the previous pin was 16.0+, itself two major versions stale).
- **Mechanism:** Subcommands, flags, and `dataformat --fields` names get renamed between major
  releases.
- **Symptom:** "Unknown flag", "field(s) [...] not recognized", or different output structure.
- **Fix:** Don't trust a version pin in any doc, this one included. Re-verify against
  `datasets <subcommand> --help` / `dataformat tsv <type> --help` for the version actually
  installed (`datasets --version`; `conda update ncbi-datasets-cli` if you need current).

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| "command not found: datasets" | Not installed | `conda install -c conda-forge ncbi-datasets-cli` |
| Subcommand not found | Old version | `conda update ncbi-datasets-cli`; re-check `datasets --help` |
| Slow 1000-genome pull | Serial download | Use `--dehydrated` + aria2c |
| "Unknown field" in dataformat | Wrong field name | Check `dataformat <type> --help` |
| Throttled bulk pull | No API key | Pass `--api-key` |
| `--reference` returns 1 per species | By design | Drop the flag or use `--assembly-level` |
| Zip checksum validation fails on `download` | Truncated or corrupt transfer | Re-run the download; persistent failure -> investigate network |
| `rehydrate` says "All N files already rehydrated" but files are wrong/tiny | It only checks that files exist, not size or checksum (see Checksum verification) | Size-check against `dataset_catalog.json`, delete mismatches, rehydrate again |
| `{"total_count": 0}`, exit code 0 | Accession doesn't exist, was withdrawn, or was superseded | Exit 0 alone is not success for `summary`/`download` -- check for a nonzero record/file count too; verify the accession at ncbi.nlm.nih.gov/datasets |
| "gene requires an at-or-below-species-level taxon" | `--taxon` given a clade (e.g. Mammalia), not a species | Use `--ortholog <clade\|all>` for cross-species gene queries instead |
| "The taxonomy name '--as-json-lines' is not exact" (unrelated taxa suggested) | Bare `--ortholog` flag swallowed the next flag as its value | Always give `--ortholog` an explicit value: `--ortholog all` or `--ortholog <taxon>` |
| `dataformat version` prints `undefined` | Unimplemented on this build | Check `dataformat --help`'s banner, or confirm it shipped with a known-good `datasets` binary |

## References

- O'Leary NA, Cox E, Holmes JB, et al. (2024) Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. *Sci Data* 11:732.
- NCBI Datasets documentation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- NCBI. Datasets CLI usage. https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/

## Related Skills

- entrez-search - For PubMed, custom queries, and non-genome data
- entrez-fetch - For single-record fetches outside genome/gene scope
- batch-downloads - Bulk E-utilities (when not genome-scale)
- sra-data - Raw sequencing reads (NOT covered by Datasets)
- ensembl-rest - Ensembl REST as alternative for Ensembl-native species
- ortholog-inference - Compara/OMA/OrthoDB for tree-aware orthology
