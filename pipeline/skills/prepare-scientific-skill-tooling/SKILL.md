---
name: prepare-scientific-skill-tooling
description: Prepare or refresh the reproducible tools, environments, public inputs, and per-surface coverage map for one normalized scientific Skill or an ecosystem batch. Use for a full pre-audit tooling pass or a focused tooling-delta pass after fixes.
---

# Prepare Scientific Skill Tooling

Make one normalized Skill executable for later workers. Do not audit its
quality, repair its instructions, or change product bytes.

Any subagent you spawn runs on Sonnet unless the brief names another model;
never Haiku.

## Load context

Verify the candidate identity, working path and branch,
lane ownership, existing `TOOLS.md`, and requested mode: `full` or `delta`.
Inspect relevant Git status, preserve all pre-existing changes, keep source
checkouts read-only, and work only on the Skills your brief assigns (one, or a batch the relay allows) and this phase only. Do not change Skill
bytes or commit, push, release, submit, or publish a product repository.

## Inventory the complete surface

Map every advertised or shipped runnable surface to its requirements:

- Python, R, shell, notebooks, compiled languages, and command-line tools;
- packages, system libraries, wrappers, runtimes, and version constraints;
- models, databases, reference files, and real public test inputs;
- remote services, authentication, licenses, registrations, hardware, and
  graphical or native-application requirements.

In a router-shaped Skill (`routes/` under the Skill root) each route is one
surface: its command, script, environment and input. Routes in the first table
of `SKILL.md` that the Skill's main tool serves are core; a route to an
alternative tool is optional.

Classify each surface as core or optional. Core surfaces are the Skill's
common end-to-end workflow and any materially distinct mode its `SKILL.md`
workflow requires; optional surfaces are alternative tools or modes it offers.
Do not build an optional surface that is heavy: it needs a GPU, a single smoke
run is expected to take over 10 minutes, its inputs exceed 5 GB, or it needs
registration or a paid license. List it in `TOOLS.md` as
`heavy-optional (not tooled)` with the reason, and record in the handoff that
the Skill must label it as not executed. Build a heavy core surface, and stage
its inputs once for the ecosystem.

Stage per ecosystem (the upstream category), not per Skill. Resolve the
ecosystem root from `F:\OpenScience\audit-envs\INDEX.md`, never by directory
name: several roots are named after an analyst (`crispr-screens` is
`crispr-screen-analyst\`). If the category has no row, create
`audit-envs\<category>\` and add the row in the same pass. Inputs sit under the
root's `public-data\`, with a README there holding one row per file: source
URL, bytes, sha256, and licence. Before any download or environment build, read
the index row, that README, and the ecosystem's existing `TOOLS.md` files,
verify a matching entry is live, and reuse it; add only what is missing. Put inputs derived from a download under
`derived\` with the script that made them. Never create a per-Skill or per-run
staging directory, and never delete staging at the end of a run.

When the brief assigns several Skills of one ecosystem, build each shared
environment and public input once and write one `TOOLS.md` per Skill that links
the shared rows.

For a delta pass, compare the fix handoff and diff with the prior inventory.
Verify unaffected surfaces retain a live environment; rebuild only what the
change invalidated. If impact is uncertain, expand the delta until every
affected surface is known.

## Build reproducibly

Use the `science` WSL distribution by default. Keep `/mnt/openscience` as its
only Windows mount and keep Windows interop disabled; never widen that boundary
for convenience. Use the existing `science` user and micromamba environments.
Stage Windows inputs under `F:\OpenScience`. Use Docker or a justified native
Windows exception when the real tool requires it. Coordinate shared mutations
with a lock, prefer isolated environments for conflicts, and snapshot relevant
versions. Kill only run-owned processes.

Prepare bounded public caches and at least one small real public dataset when
the audit would otherwise fetch them repeatedly. Keep heavy artifacts under
the disposable run root, never in the Skill package.

For each accessible primary tool, save and run a meaningful smoke test with a
bounded real public input when available. Parse data back, check scientifically
meaningful values, and render or open visual output at intended size as
appropriate. Successful installation, syntax, imports, `--help`, file
existence, or exit code zero alone is insufficient. Do not install or download
gated assets by bypassing their controls; record the exact blocker, affected
workflow, user action, and rerun steps.

## Write `TOOLS.md`

Place the Skill's `TOOLS.md` in the records repository at
`audits/skills/<skill-id>/tooling/`, not in the product Skill; shared rows link
the ecosystem `TOOLS.md` under the staging root. It
must contain:

- full or delta mode and candidate tree identity;
- one row per runnable surface (per route in a router-shaped Skill) with
  runtime/tool, environment, version, invocation, smoke evidence, input/cache
  path, and status;
- environment fingerprint and enough activation detail for a new worker;
- public data/model provenance and cache location;
- exact restricted, unavailable, or resource-infeasible blockers plus user
  action and rerun instructions;
- surfaces verified unchanged during a delta pass;
- items explicitly out of scope.

## Write the routing cases

For a router-shaped Skill, write `routing-cases.json` beside the Skill's
`TOOLS.md`: one
case per route in the first table of `SKILL.md`. The audit worker runs them
with `tools/routing_check.py`, whose header documents the format.

- `request`: one or two sentences a user would type, describing their files
  and question. Never name a route, a script, or a tool the table row does not
  name.
- `data`: a directory of small real input for that request, a few MB at most,
  cut down from the staged input and kept under the ecosystem's `derived\`
  with the script that made it. An agent that finds empty or placeholder files
  stops to report them and never reaches the route.
- `expect`: only when the route sends this request to a reference instead of a
  script; name that reference. Leave it out for a route whose command is a
  tool call shown in the route.
- `allow_before`: routes the Skill itself orders ahead of this one, such as a
  pipeline's QC route before its hit-calling route. Without it the check
  fails an agent for following the Skill's own order. The case's input must
  then pass those earlier routes: an agent that sees QC fail stops there, as
  it should, and never issues the command. An earlier route's tool must also
  be in the sandbox image, or the agent spends its steps on the missing tool.
- A route with no staged input gets one cut or simulated from real data, with
  the script that made it; say in `TOOLS.md` what is simulated.

## Hand off

Replace the canonical handoff's current state, keep it under 120 lines, and
link evidence rather than pasting output. Include `TOOLS.md`, its fingerprint,
evidence paths, blockers, and next role. Route to `audit-scientific-skill` when
no usable exact audit exists; otherwise, route to `fix-scientific-skill` or the
orchestrator's stated next phase. Do not edit the Skill, assign audit scores,
create a product commit, or continue into the next phase.
