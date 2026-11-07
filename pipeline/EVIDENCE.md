# Evidence and recovery index

Preserved on 2026-11-07. This is an archival closeout: no new scientific
certification, Skill promotion, model run, or remote publication occurred.

## Committed source and evidence

| Material | Pinned source |
| --- | --- |
| 195 published scientific Skills and provenance | [Shelf before closeout, `f162b3a9fd54e41c8f9e06bb969f86d079788518`](https://github.com/mrsonord2240/optimized-scientific-skills/tree/f162b3a9fd54e41c8f9e06bb969f86d079788518) |
| 655 audit report versions, viewers, scripts, and source records | [Records `audits/skills/`, `e9dcbbb67b1e3f02922e5d84b6778d757dce08d1`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/audits/skills) |
| Finding dispositions | [Records `fixes/`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/fixes) |
| Six pipeline Skills and shared contracts | [Original `.claude/skills/`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/.claude/skills), preserved here under [skills/](skills/) |
| 24 experimental `ca-` recuts and their actual states | [Records `skills/README.md`](https://github.com/mrsonord2240/optimizing-agent-science-skills/blob/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/skills/README.md) |
| Phase handoffs for those candidates | [Records `handoff/`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/handoff) |
| Compatibility fork | [bioSkills-Improved `9ef4d6d08244ed94b6584a0151458e09cecc727f`](https://github.com/mrsonord2240/bioSkills-Improved/tree/9ef4d6d08244ed94b6584a0151458e09cecc727f) |
| Evaluation method, all 263 result records, and recount | [EVALUATION.md](EVALUATION.md), [results.jsonl](evaluation/results.jsonl), [summary.json](evaluation/summary.json) |

Report versions are historical evidence, not 655 independently certified
Skills. Each report must be read with its `record.json`, exact source identity,
limitations, supersession information, and saved execution/assertions. The
closeout does not upgrade experimental or incomplete states.

## Local preservation archive

Directory:

```text
F:\OpenScience\archives\2026-11-07-bioskills-closeout\
```

| File | Preserves | Approximate size |
| --- | --- | ---: |
| `audit-records.bundle` | All Git refs/history of the records repository, including process tooling, audit evidence, and recuts | 31.8 MiB |
| `optimized-shelf-pre-closeout.bundle` | All Git refs/history of the published shelf before this commit | 6.0 MiB |
| `bioskills-compatibility.bundle` | All Git refs/history of the downstream fork | 20.7 MiB |
| `evaluation-evidence.zip` | 4,512 files: task source/data, three arm snapshots, original results, trajectories, events, status, interrupted attempts, and logs | 31.7 MiB |
| `workspace-notes.zip` | The fourteen environment `TOOLS.md` notes and small normalization records | See local manifest |
| `evaluation-files.json` | Per-file paths, sizes, and SHA-256 hashes for the evaluation archive | Inventory |
| `ARCHIVES.json` | Archive sizes, source commits, checksums, and verification receipts | Inventory |
| `shelf-readme-before-closeout.md` | The pre-existing README edit, including its 195-Skill count correction | Original working copy |
| `records-untracked-validate.bats` | The pre-existing untracked records test, kept outside the closeout commit | Original working copy |
| `handoff-before-closeout.md` | The former canonical handoff, retained as history | Original working copy |
| `environment-cleanup/` | Package metadata, recovered container outputs, exact removal plans, and completion receipts from the subsequent environment teardown | See [ENVIRONMENT_CLEANUP.md](ENVIRONMENT_CLEANUP.md) |

The four main archive checksums are also committed in
[local-archives.json](local-archives.json). Git verified all three bundles;
the evaluation ZIP passed CRC checks and every archived file was read back and
matched to its SHA-256. A final shelf bundle is added locally after the closeout
commit, with its receipt in the local `ARCHIVES.json`.

The ZIP excludes credentials (`.env`), Python caches, the duplicate upstream
tarball extraction, and redundant per-run `ws/` directories. It retains the
central task inputs and arm snapshots, including experimental draft bytes.
The original `F:\Bio-Evals` workspace remains intact. The archive is a local
copy on the same drive, not an off-machine backup.

## Material retained in place

The initial closeout preserved the large working trees in place. Sam then
authorized removal of the audit environments. Installed Python, Conda, and R
packages, selected scientific containers/images, and Linux package caches were
removed after their metadata and unique container outputs were preserved.
The exact scope and storage outcome are in
[ENVIRONMENT_CLEANUP.md](ENVIRONMENT_CLEANUP.md).

The large `F:\OpenScience` datasets, fix evidence, raw outputs, scripts, and
original repositories remain in place and were **not fully duplicated**. The
Git bundles preserve published reports and saved scripts, not every raw working
file. Some audit paths contain WSL links that Windows cannot resolve. The
`science` WSL distribution is retained, but its removed audit environments must
be rebuilt before rerunning their tools.

Do not delete `F:\OpenScience`, `F:\Bio-Evals`, the original repositories, or
the WSL distribution on the assumption that this compact package contains every
raw file. Machine-specific settings and unrelated working changes remain in
their original checkouts.

## Recovery

Recover a repository into a new directory with `git clone <bundle> <new-dir>`.
The bundles contain complete reachable history and do not need the existing
checkout. Extract evaluation evidence into a new directory using its saved
relative paths; verify against `evaluation-files.json` before using it.

The imported source identities are in [MANIFEST.json](MANIFEST.json). Source
and packaged hashes differ only where line endings were normalized to LF.
Result rows carry the original file hashes; the local ZIP retains exact bytes.

The active project handoff is
`C:\Users\User\Atlas-Brain\handoff\optimizing-agent-science-skills-handoff.md`.
Its closeout state supersedes the old instruction to continue recuts or wait
for Marketplace review. Historical process instructions remain available for
an explicitly reopened project.
