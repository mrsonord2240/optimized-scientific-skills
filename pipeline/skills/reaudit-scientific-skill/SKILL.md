---
name: reaudit-scientific-skill
description: Independently re-audit one fixed scientific agent Skill with risk-based executable coverage and decide candidate readiness for exact bytes. Use only for the relay's final certification phase, or its delta mode for text-only changes to certified bytes, with a fresh auditor.
---

# Re-audit Scientific Skill

Certify or reject one exact fixed candidate independently. Do not rely on the
fixer's conclusions as evidence.

Any subagent you spawn runs on Sonnet unless the brief names another model;
never Haiku.

## Establish independence and method

Verify that you did not perform the candidate's fix or initial audit. Check
origin and exact candidate identity, working path/branch, lane and audit-output
ownership, prior reports, finding dispositions, current `TOOLS.md`, and live
Git status. Preserve pre-existing changes and keep source checkouts read-only.
Work only on the Skills your brief assigns (one, or a batch the relay allows) and this phase only. Do not create a product commit, push,
release, submit, or publish.

Use the current rubric, veto rules, classification, and schema from
`skill-auditor.zip`. Treat earlier results as hypotheses to reproduce, not as
inherited passes.

Use the prepared `science` WSL environment first without widening its
`/mnt/openscience`-only, interop-disabled boundary. Use the recorded Docker or
native Windows route only where the actual tool requires it.

## Reinspect the complete tree

Repeat the full static evaluation against the exact candidate. Confirm every
prior finding's disposition, search for regressions and newly introduced
conflicts, and reconcile documentation, resources, provenance, licensing,
scientific claims, and shipped runnable surfaces.

## Execute a focused independent certification

Inventory expected runnable files, command paths, material examples, and
advertised public workflows. Independently rerun every corrected finding and
the regressions most likely to be affected, then run at least one canonical
end-to-end smoke of an untouched core workflow. Run a second end-to-end path
when it represents a materially different supported mode and is applicable.
When a program has materially different output families, ensure the combined
current evidence contains one inspected representative of each within reason.

For a router-shaped Skill, `tools/skill_preflight.py --shape` must report
`PASS`, and `tools/routing_check.py` must be rerun on the exact candidate
bytes with the lane's `routing-cases.json`; its report prints the identity it
ran against. A fixer cannot have made a case pass by editing its request:
diff the cases file against the audit's copy.

Do not rerun every prior passing input merely because it exists. Earlier saved
execution may remain part of the certification record when the auditor verifies
that the relevant Skill bytes, dependencies/runtime fingerprint, interface,
input provenance, and upstream assumptions are unchanged. Re-execute whenever
that identity is uncertain, the fix affected a shared abstraction, or earlier
evidence did not inspect meaningful output.

Inspect every result under the evidence standard. Include schemas, meaningful
values, scientific relationships, failure guards, and rendered readability as
applicable. Parsing, compilation, imports, `--help`, file existence, and exit
status never substitute for the execution and output assertions.

Do not mechanically enumerate unacceptable user inputs or duplicate upstream
validation documentation. Test an invalid or adversarial value only when silent
acceptance could plausibly corrupt scientific results, report false success,
cross a security or destructive boundary, or corrupt a public output contract.

Do not install missing tooling or download inputs; take inputs from the
ecosystem staging that `TOOLS.md` names. Route an invalidated environment or a
missing input back to a tooling-delta pass. For a paid, private, authenticated, licensed,
registration-bound, unavailable, or resource-infeasible surface, record the
exact after-action item and do not count it as executed.

## Delta mode

Use delta mode only when the brief says so and the Skill already has a
published `candidate-ready` audit. The change since those certified bytes must
be limited to prose, comments, and frontmatter, or one localized script change
of at most 20 lines. Anything larger gets the full re-audit above.

1. Verify the new identity with `tools/skill_preflight.py`.
2. Diff the candidate against the certified bytes (the optimized shelf copy or
   the certified audit's recorded path). Every change must be one the brief
   lists; report any other change as a finding.
3. Check each changed statement against the certified run's evidence or new
   evidence. Rerun the routing check when the change touched the `SKILL.md`
   route table, a route file, or a script's name or arguments. Execute any changed command or script path, including its
   intended failure guard.
4. Carry the certified report's scores forward. Re-score only the dimensions
   and assertions the change touches, and mark resolved findings.
5. Keep the certified report's schema; update identity, date, scores, and
   recommendations. Write to a new run directory and publish it as below.

A delta auditor may carry up to three Skills. It must not have fixed any of
them.

## Keep repairs exceptional

Use the same quantitative minor-repair budget as the initial-audit Skill only
for an unambiguous local defect with no dependency or scientific impact. Retain
the pre-repair observation, pin the report to resulting bytes, and rerun every
affected case plus focused independent regression. Anything larger fails the
pass and returns to a fresh fixer.

## Decide candidate readiness

Produce the complete schema-valid audit record and apply the canonical
readiness gate. Unless the current rubric is stricter, require final score at
least 85, static score at least 80, execution average at least 85, Layer 1 at
least 32, Layer 2 at least 48, assertion pass rate at least 90 percent, no veto
or open P0, and representative inspected execution for every materially
distinct required workflow or output family, using verified reusable evidence
where allowed above. A router-shaped Skill also needs every core route's
command executed, `--shape` `PASS`, and a routing check with every first-table
case `PASS`. This result makes the exact working bytes `candidate-ready`;
the orchestrator must still commit them to make them `ready` and pass intake to
make them `done`.

Publish the record and regenerate the audit index, backlog, and status outputs
at `audits/STATUS.md` and `audits/STATUS.html`. The generator derives their
counts from `audits/CORPUS.json` and published audits; never edit the counts by
hand. If it fails, resolve or report the concrete blocker before transition.

Recommendation priorities are `P0`, `P1`, or `P2` only; the records indexer
rejects anything else. A surface marked `heavy-optional` in `TOOLS.md` is
`static-only` and does not block readiness when the Skill labels it as not
executed.

For modular raw evidence, keep origin provenance and exact candidate identity
in `source-identity.json`, not as extra fields in the strict report schema.
Record the candidate's absolute `path` there, put every published script or
input under the run's `scripts/`, and use a new run directory for every run. An
uncommitted candidate needs a deterministic 64-character content SHA-256 plus
its file manifest; a 40-character Git tree id is valid only when it identifies
the exact audited subtree. Publish the explicit run root and only the named
saved scripts/inputs, then check the generated views:

```powershell
python tools/publish_audits.py --repo <records-root> --skill <skill-id> `
  --run-dir <raw-run-root> --artifact <script-or-input> [...]
npm run audits:index
npm run audits:check
```

This local records write does not authorize any push, submission, release, or
other remote mutation.

## Hand off

Replace the canonical handoff's current state, keep it under 120 lines, and
link evidence instead of pasting output. Include the certified or rejected
candidate identity, report and evidence paths, readiness metrics, every failed
or blocked surface, remaining finding IDs, worktree state, and next role. Route
failures to
`fix-scientific-skill`; route a pass to the orchestrator as `candidate-ready`.
Retire without creating product commits, running Marketplace intake, or
continuing into publication work.
