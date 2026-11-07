# Phase handoff contract

Use one canonical handoff per active Skill, also when a batched worker carries
several Skills. Replace its current-state sections;
do not append a diary. Link evidence instead of pasting it. A new worker must
be able to resume from this file, the live worktree, and linked artifacts
without chat history.

Keep the handoff under 120 lines unless a concrete blocker requires otherwise.
Use paths, hashes, finding IDs, and short facts. Do not repeat rubric text,
command output, diffs, or resolved-history narrative.

## Template

```markdown
# Handoff: <skill-id> / <next-phase>

- Updated: <ISO-8601 timestamp with zone>
- Lane: <1-5>
- Status: <ready-for-phase | blocked | phase-failed | candidate-ready>
- Owner leaving: <agent or run identity>
- Next role: <exact worker Skill name>

## Source identity

- Origin: <repository>@<full commit>:<path>
- Working tree: <absolute path>
- Branch/worktree: <name and starting commit>
- Candidate tree hash: <deterministic tree or content hash>
- Applicable audit: <path and audited identity, or none>

## Completed this phase

- <at most five outcome bullets; link evidence>

## Required next actions

1. <ordered, bounded action with finding IDs>

## Open findings and blockers

| ID | Severity | State | Evidence | Required disposition |
|---|---|---|---|---|
| <id> | <P0-P2> | <open/blocked/deferred> | <path> | <action> |

## Environment and evidence

- Tool inventory: <TOOLS.md path and fingerprint>
- Run evidence: <paths>
- Restricted-access items: <IDs or none>
- Tooling impact: <none | changed | uncertain> (<one-line reason>)

## Worktree safety

- Run-owned changes: <paths>
- Pre-existing/user-owned changes: <paths or none>
- Records state: <commit/hash or uncommitted paths>
- Product commits/pushes: none

## Transition assertion

- Next-phase prerequisites met: <yes/no>
- If no: <single precise condition>
```

Add only the receiving role's fields: structural summary and runnable surfaces
after normalization; coverage map and environment fingerprint after tooling;
report identity and finding IDs after audit; all finding dispositions and
tooling impact after fix; readiness decision and failed/blocked surfaces after
re-audit.

Reject a handoff with ambiguous candidate identity, missing evidence paths, an
unclassified tooling impact, or unexplained worktree changes.
