# Optimized Scientific Skills

Agent Skills for computational science that have been **audited against executed evidence, fixed where
the audit demonstrated a defect, and re-audited by an agent that did not write the fix**.

Nothing enters this repository until it has been through that process. A Skill that has not been
refined here lives only in its upstream repository — see [REMAINING.json](REMAINING.json).

## Attribution

These Skills originate in **[GPTomics/bioSkills](https://github.com/GPTomics/bioSkills)**, MIT
licensed, by GPTomics. That repository was archived on 2026-08-15 and accepts no issues or pull
requests, which is why refinement happens here rather than upstream.

The upstream `LICENSE` is preserved verbatim as [LICENSE](LICENSE). Every Skill remains under its
original licence. Per-Skill provenance — the upstream repository, the exact commit, the upstream path,
and whether the content was modified — is recorded in [PROVENANCE.json](PROVENANCE.json).

**We did not write these Skills.** What we added is the audit evidence, and the fixes that evidence
justified.

## What "refined" means here

| | count |
| --- | ---: |
| Skills in this repository | **75** |
| substantively modified by us, each with a fix log | 56 |
| verified unmodified apart from a declared licence | 19 |
| audit coverage | **100%** |

Every Skill here is deployable with no open P0 finding. Each carries a score and grade from a
skill-auditor run that executed the Skill's own code against real or synthetic data and recorded what
it produced — not a review of the prose.

The 19 unmodified Skills passed audit without needing a change. Their only difference from upstream is
a `license: MIT` declaration added to the frontmatter, so that the licence they already carry is
machine-readable. The audit report is what we contributed to those, not the text.

Two Skills were audited and **rejected**: `bio-experimental-design-sample-size` (67, Reject, two open
P0 findings) and `bio-experimental-design-multiple-testing` (82, not deployable). Both are listed in
`REMAINING.json` under `excluded` with their scores. They are not here, and they are not silently
omitted either.

## Layout

```
skills/<skill-id>/        SKILL.md, usually usage-guide.md, and examples/ or scripts/
PROVENANCE.json           per-Skill upstream repository, commit, path, modification status, audit result
REMAINING.json            what has not been refined yet, and what was audited and excluded
LICENSE                   the upstream MIT licence, verbatim
```

The tree is flat and the directory name equals the `SKILL.md` frontmatter `name`, so a Skill's
directory, its id, and its identity everywhere else are the same string. Upstream directory structure
is not encoded in the tree; it is recorded per Skill in `PROVENANCE.json`.

## Audit records

The full evidence — report JSON, the viewer, the scripts the auditor actually ran, and the fix log for
each modified Skill — is published separately in
[optimizing-agent-science-skills](https://github.com/mrsonord2240/optimizing-agent-science-skills)
under `audits/skills/<skill-id>/`.

## Staging

Work in progress lives in
[bioSkills-Improved](https://github.com/mrsonord2240/bioSkills-Improved), which holds the whole
upstream corpus with fixes in flight. A Skill is promoted here only once its re-audit passes.
