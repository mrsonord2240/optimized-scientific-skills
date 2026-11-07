---
name: normalize-scientific-skill
description: Normalize one scientific agent Skill, or a relay-assigned batch, before audit by resolving safe software-version drift, removing instruction redundancy, and migrating reusable code and conditional detail into routed resource directories. Use when an optimization relay assigns the normalization phase for an exact Skill.
---

# Normalize Scientific Skill

Normalize one exact Skill without performing its behavioral audit or fix pass.

Any subagent you spawn runs on Sonnet unless the brief names another model;
never Haiku.

## Start safely

Verify the assigned Skill, origin repository/path/commit, working
repository/path/branch, lane, phase, and ownership. Inspect Git status in every
repository you may touch. Treat pre-existing modified and untracked files as
user-owned; use a run-owned branch or isolated worktree and stop for any
unisolated overlap. A worktree you create is a sparse checkout holding only the
assigned Skills (`git worktree add --no-checkout`, then
`git sparse-checkout set --cone skills/<skill-id>`, then `git checkout`), never
the whole shelf. Keep provider-source checkouts read-only.

Work only on the Skills your brief assigns (one, or a batch the relay allows) and this phase only. Do not reset, clean, rebase, overwrite,
delete, or reconcile unrelated work. Do not commit a product repository, push,
open a pull request, release, submit, or publish.

## Inventory before moving

Read the worked example first, every file under its `SKILL.md`, `routes/` and
`scripts/`: `references/example/bio-differential-expression-deseq2-basics/`.
On the pseudobulk eval task that cut went from 2 of 5 to 5 of 5, with every
run reading the route and calling its script. Then read the entire tree of
the assigned Skill. Map:

- every instruction and repeated passage;
- inline and external runnable code;
- conditional, method-specific, and reference-heavy guidance;
- links among `SKILL.md`, `scripts/`, `references/`, and `assets/`;
- licenses, attribution, author information, and provider provenance.

Also compare software names, commands, APIs, formats, and version constraints
across the entire tree with the currently supported software. Treat a mismatch
caused by an ordinary software update as normalization work rather than waiting
for audit to rediscover it.

Do not infer redundancy from similar headings alone. Preserve distinct
conditions, caveats, and scientific meaning.

## Normalize structure

- Cut the Skill to a router that gets more detailed as it branches. Fidelity
  to the provider's text and layout is not a goal (Sam, 2026-10-04); correct
  science, credit and licence are.
  - `SKILL.md`: a table from what the user has to one route file, the few
    rules that hold on every route, and nothing else. Tell the agent to read
    one route and run its command before writing code of its own.
  - `routes/<task>.md`: the command first, then only the rules that change
    what the agent does on that route, then a "Done when" line.
  - `scripts/`: one runnable entry point per route, which takes the user's
    files and writes one results file. Build traps into the script (refuse a
    missing reference level, stop on an un-normalized matrix) instead of
    warning about them in prose. A script that serves one route takes
    the route's name (`routes/pseudobulk.md` runs `scripts/pseudobulk_de.R`);
    one that serves several is named for what it runs (`deseq2_de.R`). An
    eval agent that skipped the route file guessed a script name from the
    route name, found nothing, and invented its output.
  - `references/`: per-tool and per-topic detail a route points to by name.
  - Work that belongs to a sibling Skill gets one route: a table from the
    design or need to that Skill's name, and nothing that repeats it.
- Start with the point. The first thing under the title is the route table or
  the default command, within 15 lines of the title. Tested versions shrink to
  one line at the end of the file.
- `SKILL.md` and route files hold only what an agent needs to do the task:
  what to run, what to choose between, when to stop, and what breaks. Move the
  rest to `references/background.md`: essays on why a method is right, tool
  taxonomies that describe mechanisms, history, and the provenance of a number
  ("measured on one 4 v 4 set"). Nothing routes there for doing the task.
- Keep a reason only when it changes what the agent does, and keep it to one
  clause: "sum raw counts, never normalized values" stays; a paragraph on why
  pseudoreplication inflates false discoveries goes. Test each sentence by
  deleting it: if the agent would act the same, it moves. Sam, 2026-10-04:
  too much language justifying techniques and giving provenance, so a reader
  reaches line 57 before the Skill gets to the point.
- Move substantial reusable executable code to `scripts/`.
- Move conditional, method-specific, or large lookup material to
  `references/` when it has a real routing purpose.
- Keep source assets in `assets/` only when they belong in generated output.
- Remove duplicate material after ensuring one authoritative copy remains.
- Repair relative links and references made stale by the moves.
- Preserve the frontmatter name, license obligations, attribution and
  provenance. Keep every analysis the Skill could do reachable from a route;
  drop human-facing guides and demos on simulated data that a script replaces.
- Resolve version discrepancies throughout `SKILL.md`, routed references,
  scripts, examples, and dependency metadata. Default to the newer supported
  software version when live documentation or tooling confirms it and the
  update causes no evident breakage, dependency conflict, format incompatibility,
  or change in scientific meaning. Update connected commands and claims
  together. Record a genuinely ambiguous or breaking migration for audit or
  user decision instead of guessing.

- Make the frontmatter Marketplace-complete: keep `name`, `description`, and
  `license`; add exactly one `category:` from `Academic Writing`,
  `Data Analysis`, `Evidence Insight`, `Other`, `Protocol Design`, chosen from
  what the Skill does; add `author:` with the original author credit from the
  provider (for GPTomics bioSkills: `author: GPTomics`).
- Reduce the frontmatter `description` to its trigger: one sentence of the
  form `Use when <task or situation>.` Remove tool and method lists, coverage
  summaries, and "for X see <other Skill>" pointers; keep that material in the
  body. Name the domain and task precisely enough to tell the Skill apart from
  its siblings in the same category, and state nothing the Skill does not do.
  In a small trial (Sam, 2026-10-03) trimmed descriptions made agents more
  likely to select a Skill, select the right one, and finish the task.
- Move literature citations out of `SKILL.md` and every routed instruction
  file into `references/citations.md`: the `## References` section, and inline
  author-year or journal citations in prose and table cells. Keep the
  instruction, threshold, or number the citation supported and delete only the
  attribution; never drop a rule because its source moved. List each source
  once in `citations.md` beside the claim it supports. Nothing routes to that
  file for doing the task; it is evidence for reviewers. A method's own name
  (Benjamini-Hochberg, Wald test) is not a citation and stays. Sam,
  2026-10-04: citations do nothing toward getting the task done and crowd the
  text an agent has to read.
- Write every text file as UTF-8 without BOM, with LF line endings. Remove
  `__pycache__`, `*.pyc`, dot-prefixed paths, and LICENSE copies below the
  Skill root.

Do not redesign scientific methods, change analytical defaults, add analyses
the Skill did not offer, install tooling, or resolve behavioral findings
unrelated to safe version normalization in this phase. Record such issues for
audit. A route's entry-point script that wraps the Skill's existing method is
structure, not a new analysis; list every script you wrote in the handoff.
Do not run them: tooling smoke-tests every route and audit runs every core one.

## Verify the normalized tree

Check that every analysis the Skill offered is reachable from the route table.
Leave behavioral audit of the science to tooling and audit workers.

Then run, from the records repository root,
`python tools/skill_preflight.py --shape <normalized-skill-dir>`. It must
report `PASS`. `--shape` fails a tree with no route table within 15 lines of
the title, a route nothing names, a named `routes/`, `scripts/` or
`references/` path that does not exist, a citation left in `SKILL.md` or a
route, or a description that is not a `Use when` trigger. Take the candidate
identity from its output; never compute or copy an identity by hand. Carry any
`warn` lines into the handoff.

Summarize migration at file or purpose level. Do not catalog every moved or
deleted passage; the diff is that record.

## Hand off

Replace the canonical handoff's current-state sections; do not append a diary
or paste evidence. Keep it under 120 lines and link durable paths. Include the
resulting candidate identity, concise structural summary, the scripts you
wrote, open ambiguities, and all run-owned versus pre-existing changes. Set the next role to
`prepare-scientific-skill-tooling`. Do not create a product commit, push, or
continue into tooling.
