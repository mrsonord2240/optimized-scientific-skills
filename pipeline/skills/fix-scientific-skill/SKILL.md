---
name: fix-scientific-skill
description: Repair one audited scientific agent Skill against a durable finding ledger, execute changed behavior, and classify tooling impact for independent re-audit. Use for the relay's potentially long, resumable fix phase, or a text-only fix batch.
---

# Fix Scientific Skill

Resolve the assigned audit findings without scoring or certifying your own
work. Keep the phase resumable even when it spans many turns.

Any subagent you spawn runs on Sonnet unless the brief names another model;
never Haiku.

## Load the exact assignment

Read the complete active audit report, ordered finding ledger, `TOOLS.md`, and
live diff. Verify the origin and candidate identities, working path/branch,
lane ownership, and audit identity. Inspect relevant Git status, preserve all
pre-existing changes, keep source checkouts read-only, and isolate unexplained
changes before editing. Work only on the Skills your brief assigns (one, or a batch the relay allows) and this phase only. Do not create a
product commit, push, release, submit, or publish.

Do not repeat normalization as a ritual. Preserve its resource routing and do
not reintroduce duplication. When a finding genuinely requires structure to
change, keep `SKILL.md` routing correct and summarize the result; the diff
remains the detailed migration record.

## Maintain one finding ledger

Use the audit's stable IDs. For every inherited finding, maintain exactly one
state: `open`, `in-progress`, `fixed`, `not-reproduced`, `blocked`, or
`deferred-with-rationale`. Link the implementing diff and verification evidence
instead of copying their contents into prose. Add newly discovered defects with
new IDs and severity.

Fix every known P0 by default. Fix a P1 or P2 only when it blocks readiness,
produces wrong or unusable output, invalidates an advertised workflow, or is a
text-only correction; record every other one as `deferred-with-rationale`
("not needed for readiness"). Do not work toward a higher score than the
readiness gate. Surface wider product, scientific, license, access, or
architectural decisions instead of guessing.

A heavy optional surface the tooling worker left untooled gets a short label in
the Skill saying it was not executed here; do not tool or execute it.

### Text-only batch mode

When the brief assigns a batch of `candidate-ready` Skills with text-only P2s,
change only prose, comments, and frontmatter; no runnable behavior may change.
For each Skill, record the finding IDs, changed files, and new identity in its
handoff and fix log, and route it to `reaudit-scientific-skill` in delta mode.
If a finding turns out to need a runnable change, leave it open and say so.

In a router-shaped Skill a finding belongs to a route. Fix a routing-check
finding in the Skill: the table row, the route file, or a missing entry point.
Never edit a case's request to lead the agent there.

## Implement deliberately

Prefer the smallest durable correction that addresses the underlying cause.
Preserve provenance, licensing, trigger accuracy, supported interfaces, and
scientific intent. Do not mask failures, weaken assertions, narrow advertised
scope without evidence, fabricate outputs, or substitute a materially
different workflow merely to obtain a pass.

Whenever a claim is corrected, use `rg` across the complete Skill tree for the
old wording, close variants, duplicated examples, and connected statements
whose validity depended on it. Correct the whole invalidated claim family in
the same pass or record each intentional survivor with evidence. Do not stop at
the line named by the finding.

Use the prepared tooling and take inputs from the ecosystem staging that
`TOOLS.md` names; do not download into the run directory. If the fix needs a new dependency, runtime, version,
wrapper, model, dataset, service, or executable path, record that fact and do
only the minimum safe implementation work needed; environment installation and
validation belongs to a subsequent tooling-delta worker.

## Verify changed behavior

Run syntax or import checks as preflight, then execute every changed accessible
behavior with meaningful inputs. Save commands and versions; parse structured
outputs, check meaningful scientific values and relationships, and inspect
rendered readability where applicable. Exit code zero or file existence alone
is insufficient. Rerun focused regressions for adjacent behavior affected by
the change. A change to a script or reference that several routes share reruns
each of those routes' commands. Reuse immutable evidence for an unchanged surface when its bytes,
dependency/runtime fingerprint, interface, and relevant upstream assumptions
still match; do not replay expensive unrelated workflows ceremonially. A fix
may remain blocked when access is unavailable, but it must not be marked
verified. Do not bypass authentication, payment, licenses, or registration.

Long phases must update the same finding ledger and canonical handoff after
each coherent milestone. Keep evidence in files and replace current-state
sections rather than appending narrative. Records-repository checkpoint commits
are allowed when they contain only run-owned resumable state; product commits
are not.

## Classify tooling impact

Before handoff, set one value:

- `none`: no dependency, runtime, version, wrapper, model, dataset, service,
  executable path, or runnable surface changed;
- `changed`: at least one of those changed, with affected surface IDs;
- `uncertain`: evidence is insufficient to prove impact is absent.

Do not use `none` merely because the old environment still starts.

## Hand off

Run `python tools/skill_preflight.py --offline <skill-dir>` from the records
repository root; it must report `PASS`, and its identity is the candidate
identity for the handoff.

Replace the canonical handoff's current state, keep it under 120 lines, and
link evidence rather than pasting output. Include every finding disposition,
changed file, execution record, remaining blocker, exact candidate identity,
worktree safety, and tooling impact. Route `changed` or `uncertain` to
`prepare-scientific-skill-tooling` in delta mode; route `none` to
`reaudit-scientific-skill`. Retire without assigning a readiness score,
creating a product commit, or continuing into final audit.
