---
name: audit-scientific-skill
description: Conduct a bounded diagnostic audit of one normalized scientific agent Skill (or a relay-assigned batch of light Skills) using prepared tooling, producing schema-valid findings and executable evidence without claiming final certification. Use for the relay's initial audit phase when no usable exact audit exists.
---

# Audit Scientific Skill

Audit one exact normalized candidate to identify and prioritize work. This is a
comprehensive but bounded diagnostic pass. Any subagent you spawn runs on
Sonnet unless the brief names another model; never Haiku. Independent re-audit later focuses
on corrections, affected regressions, and representative core workflows.

## Load context and method

Verify the origin and exact candidate identity, working path and branch, lane
and audit-output ownership, current `TOOLS.md`, and separation from other
workers.
Inspect relevant Git status, preserve all pre-existing changes, and keep source
checkouts read-only. Work only on the Skills your brief assigns (one, or a batch the relay allows) and this phase only. Do not make a
product commit, push, open a pull request, release, submit, or publish.

Use the current `skill-auditor.zip` rubric, veto rules, classification, and
report schema from the records repository. Extract it to the disposable run
area if needed. Never silently substitute a remembered rubric.

Use the prepared `science` WSL environment first without widening its
`/mnt/openscience`-only, interop-disabled boundary. Use recorded Docker or
native Windows exceptions only when required by the actual tool.

## Inspect statically

Read the complete Skill tree and evaluate frontmatter, trigger precision,
instruction consistency, progressive disclosure, resource routing,
provenance, licenses, safety boundaries, scientific-method claims, examples,
and correspondence between documented and shipped runnable surfaces.

Treat normalization as already complete. Report missed structural problems;
do not restart a broad migration inside the audit.

## Execute representative end-to-end workflows

Use the prepared environment without installing a new stack, and take inputs
from the ecosystem staging that `TOOLS.md` names. Do not download into the run
directory; record a missing input in the handoff for a tooling-delta pass. Run the most
common public workflow end to end, from realistic input through inspected
scientific output. Run a second end-to-end workflow when it exercises a
materially different supported mode, runtime, algorithm, or integration and a
bounded case is available. When a program produces multiple materially
different output families, generate and inspect one representative of each
within reason. Execute any migrated executable whose preservation is uncertain.

For a router-shaped Skill, run every core route's command as its route file
shows it, on the staged input, and trip each trap a script builds in once (a
missing reference level, an un-normalized matrix). The stop message must name
the fix. An optional route follows the rules for optional surfaces.

Save commands, inputs, versions, outputs, and assertions. Parse structured
output and meaningful values; inspect scientific relationships and rendered
readability where applicable. Syntax, imports, file existence, and exit code
zero are preflight only.

Do not mechanically enumerate unacceptable user inputs or reproduce validation
rules already maintained in the upstream software's documentation. Exercise an
adversarial or invalid input only when silent acceptance could plausibly corrupt
a scientific conclusion, report false success, cross a security or destructive
boundary, or corrupt a public output contract. Ordinary input mistakes are not
a required test matrix.

A surface marked `heavy-optional (not tooled)` in `TOOLS.md` is `static-only`
with reason `heavy-optional`; check only that the Skill labels it as not
executed, and file a finding when it does not.

You may defer expensive secondary examples and redundant permutations. Mark a
materially distinct advertised workflow that lacks representative evidence as
`static-only` or `blocked`; do not treat a combinatorial parameter or invalid-
input matrix as a separate workflow. Do not score documentation or source
inspection as execution.

If missing tooling prevents meaningful diagnosis, stop that surface and route
the lane to `prepare-scientific-skill-tooling`; do not install dependencies
during scoring.

## Run the routing check

For a router-shaped Skill, run from the records repository root:

```powershell
python tools/routing_check.py <candidate-skill-dir> <routing-cases.json> --out <run-root>\routing
```

A cheap model gets each case's request and input and the Skill; the check
passes a case when the model opens the right route first and then runs its
script. It executes nothing from the Skill, costs a few cents, and needs
Docker; when containers will not start, run `docker desktop restart` and
rerun (Sam's standing permission).

- `FAIL`: read the failing runs under `--out` and file a P1 finding on that
  route that says why the agent went elsewhere: the table row does not
  describe the request, the route needs a file the user does not have, or the
  command is not first.
- `ERROR`: infrastructure; rerun those cases with `--routes a.md,b.md` before
  judging. HTTP 402 means the OpenRouter account is out of credit: rerun the
  whole check with `--model cli:claude-haiku-4-5-20251001`, which uses the
  local Claude subscription, and record the model in the report.
- A route that fails only because the agent followed the Skill's own stated
  order (QC before hit calling) is a case defect: it needs `allow_before`
  from the tooling worker, not a finding on the Skill.
- A case whose request leads the agent, or whose input does not fit the
  route, goes back to the tooling worker.

Keep `routing.json` with the run evidence.

## Use the minor-repair budget narrowly

An audit-local repair is allowed only when all are true:

- the intended correction is unambiguous and changes no scientific method,
  analytical default, claim, interface, security boundary, or dependency;
- it is one localized defect affecting no more than two files and 20 changed
  lines, excluding purely mechanical formatting;
- it requires no architectural, product, licensing, or access judgment;
- every affected case and a focused regression can be rerun in this pass.

Examples include an unmistakable typo, missing import, stale documented flag
with one authoritative replacement, or incorrect local relative path. If any
condition is doubtful, file a finding for the fixer. Pin the audit to the
post-repair bytes and retain the pre-repair observation and rerun evidence.

## Produce and publish the record

Write schema-valid `report.json`, `viewer.md`, saved run scripts and inputs,
source identity, execution classifications, findings with severity and stable
IDs, minor-repair evidence, and restricted-access after-action items. Keep a
single ordered finding ledger for the fixer. Recommendation priorities are
`P0`, `P1`, or `P2` only; the records indexer rejects anything else.

For a modular run root, `source-identity.json` must keep the immutable origin
repository/commit/path separate from the exact working candidate. Identify an
uncommitted candidate with a deterministic 64-character content SHA-256 (and
its file manifest); use a 40-character Git tree id only when it identifies the
exact audited subtree. Take that identity from `tools/skill_preflight.py`, and
record the candidate's absolute `path` in `source-identity.json`. Put every
published script or input under the run's `scripts/`, and use a new run
directory for every run because published versions are immutable. Never put
publication metadata into the strict audit report merely to satisfy older
tooling.

Publish through the records repository tooling and regenerate the audit index,
backlog, and status outputs at `audits/STATUS.md` and `audits/STATUS.html`. The
generator derives their counts from `audits/CORPUS.json` and published audit
records; never edit the counts by hand. A diagnostic score does not make the
Skill `candidate-ready` or `ready`.

Use the explicit modular mode and name every saved script or bounded input that
belongs in the public record:

```powershell
python tools/publish_audits.py --repo <records-root> --skill <skill-id> `
  --run-dir <raw-run-root> --artifact <script-or-input> [...]
npm run audits:index
npm run audits:check
```

Local records publication is required audit evidence. It is not authorization
to push, submit, release, or otherwise mutate remote state.

## Hand off

Replace the canonical handoff's current state, keep it under 120 lines, and
link rather than paste evidence. Include the report path and exact identity,
ordered open finding IDs, deferred and blocked surfaces, environment
fingerprint, repair details, worktree state, and next role. Normally route
findings to `fix-scientific-skill`; a clean pass routes to
`reaudit-scientific-skill`. Retire without creating a product commit or
continuing into fixes.
