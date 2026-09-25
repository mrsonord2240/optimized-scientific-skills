# Preliminary Marketplace submission: 20-Skill pilot

This directory contains the provider-side material for the first 20-Skill
Marketplace pilot. The source packages are pinned to commit
`524cbb7ce83008e2e037bed0bb1cd04fc61dcb82` of
`https://github.com/mrsonord2240/optimized-scientific-skills`.

The files have deliberately separate responsibilities:

- Each Skill directory contains the strict Marketplace `release.config.json`.
  Do not add provider-only fields to those files; the Marketplace schema rejects
  unknown fields.
- `intake-review-input.json` is the unapproved output of the Marketplace intake
  command. It records the exact manifest/content hashes, file counts, package
  sizes, license expression, and license-evidence hashes for the pinned source
  commit. It is evidence input, not approval.
- `pilot-review-evidence.json` records provenance, reproducible external-audit
  links, runtime dependencies, declared input/output contracts, bundled
  executables, review scope, unresolved limitations, and provider check results.
  It also retains the verified intake metrics per Skill under `package_evidence`,
  so those facts do not depend on a local generated intake file being committed.

## Attribution and provenance

The original source is GPTomics `bioSkills` at commit
`d91ed3d563019e649dc854c56ccd62551359488a`. GPTomics remains the author named
in every pilot Skill. Every pilot Skill is marked `modified`; the modifications
made in the `mrsonord2240` repositories are not represented as AIPOCH-original
content. Per-Skill upstream paths and modification evidence are recorded in
`PROVENANCE.json` and linked from `pilot-review-evidence.json`.

The external audit repository is pinned to
`https://github.com/mrsonord2240/optimizing-agent-science-skills` commit
`0d1d6c1f2e8c3ec89316b58269a4a8c795690cf5`. Its scores are retained only as
upstream audit evidence. They are not AIPOCH review scores. Each Skill now links
the exact matching report and record, with their SHA-256 hashes and audited
source identity. Those reports refer to commits in `mrsonord2240/bioSkills`;
byte equivalence with the submitted packages has not been established. Sixteen
reports explicitly say `auditor_independent: false`; four do not specify it.

## Review gate

Every `aipoch_review` record is intentionally pending, with `reviewed_by`,
`reviewed_on`, and `score` set to `null`. Those values must be filled only by a
Marketplace maintainer after reviewing the exact pinned bytes, license scope,
attribution, dependencies, and semantic contract. The official reviewed batch
JSON should be generated from `intake-review-input.json`; this provider evidence
file is not a substitute for that Marketplace review record.

Separate `provider_review` records identify Codex as the provider-side reviewer,
date the review, and state its completed scope and execution coverage. These
records do not claim to be a Marketplace maintainer's approval. The review used
separate spec and standards agents to check the prepared evidence.

## Current requirement status

1. Release configs: complete for all 20 pilot Skills.
2. Required release fields: complete and accepted by Marketplace intake.
3. Package and license hashes: complete in `intake-review-input.json` for the
   pinned source commit.
4. Review records: provider reviewer, date, scope, results and limitations are
   recorded per Skill; AIPOCH reviewer/date/score remain pending.
5. GPTomics provenance and modification attribution: recorded.
6. External audit repository: repository, commit, audit path, and fix-log path
   are pinned per Skill.
7. Runtime dependencies: cross-checked against bundled scripts and documented
   routes, including Python/R libraries, Bash, Java/Perl, command-line tools,
   reference resources, downloads and models. This is not a tested environment
   lock or evidence of clean installation of every optional route.
8. Input/output/limitation/script cross-check: static checks completed for all
   20; selected synthetic execution covers 13. Three Skills have open source or
   runtime findings below. Full functional verification remains incomplete.

## Local syntax checks

The bundled executable inventory passed Python AST parsing (52 files), Bash
syntax checking (6 files), and R parsing (10 files). The R check used
`C:\R\bin\x64\Rscript.exe` with R 4.6.1 and `--vanilla` on 2026-09-24.
The earlier PATH-only lookup did not find Rscript; that was a discovery failure,
not a missing installation. R emitted locale startup warnings for `C.UTF-8`,
but all ten parse commands exited successfully.

The verifier confirms every local packaged file equals its pinned Git blob,
allowing these syntax results to be associated with the submitted source.
Syntax checks alone do not establish functional correctness.

## Functional checks and open findings

Sixteen offline checks across 13 Skills produced 15 passes and one failure.
Each run used an isolated copy of its package. Full command arguments, source
script hash, stdout, stderr and exit code are stored under each Skill's
`provider_review.runtime_checks`. Python demos, Scanpy preprocessing/clustering,
R MaxLFQ, BAM filtering/validation, VCF normalization and germline hard-filtering
passed their stated checks. This exercises only the named routes and fixtures.

The remaining findings are explicit holds for review:

- `bio-proteomics-data-import`: QFeatures produced the expected cleaned matrix
  summary, then Windows R 4.6.1 exited with `3221225477` (`0xC0000005`, access
  violation). WSL R is available but does not have QFeatures installed. The
  Python route passed. The source of the native crash remains undiagnosed.
- `bio-single-cell-cell-annotation`: the bundled CellTypist example does not
  validate its CP10K-log1p/UMAP prerequisites and does not use the documented
  seeded over-clustering setup. It also needs an external model. This route
  was not executed.
- `bio-vcf-statistics`: the Python example skips valid `QUAL=0` values when
  computing the mean and counts only the first ALT for multiallelic Ti/Tv.
  These are source findings; the example was not executed because cyvcf2 is
  absent from the checked environments.

The variant-filtering shell example's SNV/indel-only boundary is also recorded:
it warns that MNP/mixed records are dropped. Its scoped synthetic test passed.

No source Skill files were changed in this review, so the submission remains
pinned to the original selected bytes. Resolving a source finding requires a
new source commit, refreshed release configs/intake hashes and new review.

## Reproduce provider checks

Run from this provider repository with trusted local clones available:

```powershell
node authoring/verify-pilot.mjs --marketplace F:/OpenScience/marketplace-intake/openscience-skill-marketplace --audit F:/optimizing-agent-science-skills
node authoring/smoke-pilot.mjs --python python --rscript C:/R/bin/x64/Rscript.exe
```

The verifier re-runs official intake, compares committed `package_evidence`,
optionally compares `intake-review-input.json` when it is present locally,
checks every executable and package byte, and verifies pinned audit links. It
prints JSON without executing Skill code. The smoke runner creates unique
temporary package copies, keeps its logs, and exits nonzero on any failed check.
Use repeated `--only <skill-id>` options to restrict a repeat run.

For the four Linux route tests, export packages to a drive mounted in WSL:

```powershell
node authoring/smoke-pilot.mjs --export-only --scratch-parent F:/ --only bio-alignment-filtering --only bio-alignment-validation --only bio-variant-normalization --only bio-variant-calling-filtering-best-practices
```

Pass the printed scratch directory (converted to its `/mnt/f/...` path) to
`python3 -B /mnt/f/optimized-scientific-skills/authoring/smoke-pilot-linux.py`
inside the `science` WSL distribution. That runner verifies the exported file
hashes before executing its synthetic checks. Requires pysam, samtools,
bcftools and Bash. No external datasets or services are used.

## Marketplace validation scope

At Marketplace tooling commit `e751844be4ade3e4d997ea4dbc63a0b81cf6ff75`,
`intake:skill` successfully regenerated all 20 records byte-for-byte at the JSON
data level. `validate` and all 98 tests passed. `publish:dry-run` passed its
built-in fixture rehearsal. The latter does **not** publish, approve, enroll or
build this pilot: an approved maintainer review is still required for that.
