# Audit environment cleanup - 2026-11-07

After the preservation closeout, Sam authorized removing the audit
environments. The selected installations and package caches have been removed.
No scientific audits or model evaluations were rerun.

## Removed

| Surface | Completed scope |
| --- | --- |
| Windows, under `F:\OpenScience\audit-envs` | Installed components of 46 Python environments, 19 R package libraries, and one staged Conda environment; 235 explicit targets |
| `science` WSL | Installed components of 114 scientific Conda environments and recognized packages in `/home/sci/R-lib`; 1,529 explicit targets |
| Linux package caches | 4,518 extracted packages, 4,518 package downloads, and two channel index caches under the Micromamba/Miniconda package directories |
| Docker containers | The stopped `ccr2` and `ccr3` containers |
| Docker images | `bioeval-env:latest`, `vaca-audit:1`, `pinellolab/crispresso2:latest`, and `bioconductor/bioconductor_docker:RELEASE_3_21` |

Removal receipts report no failures. Native protected paths and Linux retained
extras were checked after deletion. Package-owned components were selected
explicitly because environment directories also contained evidence and data.
Empty environment roots were removed only when empty.

## Preserved

- The 195 scientific Skills, six pipeline Skills, evaluation records, audit
  reports, 24 experimental recuts, original repositories, and closeout archives.
- Datasets, scripts, logs, raw outputs, fix evidence, and the original
  `F:\Bio-Evals` workspace. Remaining `TOOLS.md` files are historical recipes;
  their former installed environments are no longer available.
- Package identities and installation records: Python `pyvenv.cfg`, pip
  metadata/records, R descriptions, Conda metadata/history, and Docker image
  identities/build history. These support reconstruction; they are not binary
  snapshots or a guarantee that upstream packages remain obtainable.
- Unique `ccr3` outputs, inputs, R scripts, and package description, copied out
  before container removal. Both removed containers' logs are retained.
- The WSL distribution, its OS/base installation, the general
  `pdf-render-poppler` environment, `/home/sci` evidence/configuration, and
  unclassified installation leftovers. The distribution is stopped.
- Other Docker services and images. GBrain services remained running.

No whole audit category, data tree, WSL distribution, or Docker volume was
deleted. The retained raw working files are not all backed up by the compact
archives. Small unclassified leftovers and standalone tools remain for any
later, separately scoped cleanup.

## Evidence location and verification

```text
F:\OpenScience\archives\2026-11-07-bioskills-closeout\environment-cleanup\
```

| Archive | Files | Contents |
| --- | ---: | --- |
| `windows-environment-metadata.zip` | 9,442 | Native package records, environment markers, and setup notes |
| `linux-environment-metadata.zip` | 42,766 | Conda/pip/R records for the removed Linux installations |
| `ccr3-evidence.zip` | 15 | Unique container inputs, outputs, and source |

Every archive member was checked against its saved SHA-256 manifest. The local
directory also contains inventories, process checks, exact removal plans,
per-target journals, completion receipts, and the compaction result. Archive
checksums and receipt identities are committed in
[environment-cleanup.json](environment-cleanup.json); the full local archive
index is `..\ARCHIVES.json`. Original closeout archive hashes were rechecked
after cleanup. These copies are on the same drive, not off-machine backups.

## Storage outcome

F: free space increased by **69.08 GiB during native environment removal**.
Subsequent metadata preservation used some of that space; this is an observed
phase delta, not an exact final net saving or a Docker storage estimate.

Linux environment removal freed **118.56 GiB** inside WSL, and package cache
removal freed another **85.34 GiB**: approximately **203.89 GiB** in total.
Filesystem trim succeeded. The distribution was stopped before attempting
compaction, but `Optimize-VHD` was denied because the Windows session lacked
the required administrator permission. The virtual disk remains
275,840,499,712 bytes (256.9 GiB), so the internal WSL saving has **not yet been
returned to F:**. Compaction from an appropriately elevated session remains
the next storage step; it does not require deleting the distribution.

The cleanup is documented locally and does not publish or push the repository.
