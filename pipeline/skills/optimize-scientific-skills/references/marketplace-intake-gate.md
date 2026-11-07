# Marketplace intake gate

Run this local compatibility gate after assembling the unpublished optimized
batch commit and before calling any Skill `done`. It supplements, and never
replaces, the scientific re-audit.

## Boundary

Use the Marketplace repository at
`F:\optimizing-agent-science-skills\marketplace\intake\openscience-skill-marketplace`
and its `npm run --silent intake:skill` command, at the current upstream
validator: run `git fetch origin`, then run intake from a detached worktree of
`origin/main` (`git worktree add --detach F:\OpenScience\marketplace-validator
origin/main`, then `npm ci` there; move an existing one with `git -C
<worktree> checkout --detach origin/main`). Record that commit in the receipt. Intake reads committed Git
blobs from a full commit SHA and ignores working-tree changes; the manifest
must therefore pin the exact unpublished optimized batch commit.

Do not supply review or output-build arguments, add reviewer identity, build a
catalog or bundle, copy packages into a submission queue, create submission
records, register a release, enroll, publish, or push.

## Manifest inputs

In a disposable directory under `F:\OpenScience`, create one temporary
`release.config.json` per Skill. Derive rather than guess:

- ID and category from `SKILL.md` frontmatter (`name`, `category`; the
  preflight already enforces one of the five Marketplace labels);
- version from release metadata, `1.0.0` for a first release;
- canonical optimized repository URL;
- full optimized commit SHA and `skills/<skill-id>` source path;
- actual nonempty license files at that commit.

An unresolved identity, version, category, or license judgment is a user-action
blocker.

## Run and verify

From the Marketplace repository root, pass matching manifest and source
arguments for every candidate. The source is
`F:/optimized-scientific-skills`. Require:

- exit code zero;
- stdout that parses as a nonempty JSON map;
- exactly one `<id>@<version>` entry per candidate;
- expected repository, full commit, and path in every entry;
- no hidden rejection or partial failure on stderr.

Save a compact receipt in the records repository: validator commit, command
form, optimized commit, ID/version, manifest hash, returned evidence/content
hashes, timestamp, and result. Temporary review-input JSON is not approval.

## Route failures

- Manifest-only: correct the temporary manifest and rerun; unchanged Skill
  bytes do not require re-audit.
- Authoritative provider metadata or license evidence: correct it, recreate the
  unpublished commit if committed metadata changed, and rerun.
- Skill package bytes or structure: return to fix, fresh executable re-audit,
  commit recreation or amendment, and intake rerun.
- Unresolved access, license, identity, or policy judgment: exclude the Skill,
  preserve its state, and report the user action required.

Do not preserve superseded provisional commits in final product history.
