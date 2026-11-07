# Repository and worktree contract

## Repository roles

- `F:\optimized-scientific-skills` is the canonical editing source and final
  cross-source product shelf.
- `F:\optimizing-agent-science-skills` holds control state, audit and fix
  evidence, generated status views, and resumable work records.
- `F:\OpenScience` holds disposable worktrees, environments, cached public
  inputs, run outputs, and other reproducible heavy artifacts.
- `mrsonord2240/bioSkills-Improved` is a downstream, bioSkills-shaped
  compatibility fork generated one way from completed canonical bytes by
  `tools/export_bioskills.py`. Its checkout is `F:\OpenScience\bioSkills-Improved`.
- Original provider repositories are provenance sources. Keep their source
  checkouts read-only, including generated caches and interpreter artifacts.

Use user-supplied roots when they explicitly replace these defaults.

## Identify and claim the work

Before mutation, record the canonical Skill ID, origin repository/path/commit,
working path/branch/starting commit, applicable audit identity, current phase,
owner, and durable handoff path. The directory name and frontmatter `name` must
agree. Claim one phase only, for one Skill or for the batch the relay batching
table allows.

## Protect existing work

Inspect Git status in every repository the phase might touch. Treat all
pre-existing modified and untracked files as user-owned. Use a run-owned topic
branch or isolated worktree for mutation. Do not reset, clean, rebase,
overwrite, delete, or reconcile unrelated work. Stop and report an overlap
that cannot be isolated safely.

Make only changes needed for the assigned phase. Never bypass authentication,
payment, licensing, registration, or other access controls.

## Keep worktrees sparse

A worktree holds only the Skills being optimized in it, never a copy of the
whole shelf. Create it without a checkout and set the cone before populating:

```powershell
git -C F:\optimized-scientific-skills worktree add --no-checkout -b <branch> F:\OpenScience\wt\<name> <start-commit>
git -C F:\OpenScience\wt\<name> sparse-checkout set --cone skills/<skill-id> [...]
git -C F:\OpenScience\wt\<name> checkout
```

Cone mode keeps the repository's root files and the named Skill directories
only. Add a second Skill to the cone only when the same optimization needs to
edit or read it. A worker that finds a full checkout reports it; it does not
copy other Skills in.

## Product and publication boundaries

Phase workers do not commit product repositories. The orchestrator forms at
most one completed-Skill batch commit per affected product repository when the
run closes. A one-Skill batch is valid; an empty commit is not.

Workers may save run-owned evidence and resumable state in the records
repository. They must not push, open a pull request, release, submit, enroll,
publish, modify remote state, or perform Marketplace review or bundle building
without separate authorization.

Finished trees retain exact-byte identity from final audit evidence to the
optimized shelf. A post-audit Skill-byte change returns to re-audit. The
maintained bioSkills fork receives only the accepted bio-derived subset and is
never an alternate editing source.
