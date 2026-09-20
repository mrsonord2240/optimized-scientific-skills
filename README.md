# Optimized Scientific Skills

Agent Skills for computational science that have been **audited against executed evidence**. Where the
audit demonstrated a defect, the Skill is fixed and then re-audited by an agent that did not write the fix.

Nothing enters this repository without an audit that found it deployable with no open P0 finding. **Not
every Skill here has finished the process yet.** Each one's status is in [PROVENANCE.json](PROVENANCE.json)
(`fix_pass`, `reaudit`), and both lists are at the top of [REMAINING.md](REMAINING.md). A Skill that has not
been audited lives only in its upstream repository.

## Attribution

These Skills originate in **[GPTomics/bioSkills](https://github.com/GPTomics/bioSkills)**, MIT
licensed, by GPTomics. That repository was archived on 2026-08-15 and accepts no issues or pull
requests, which is why refinement happens here rather than upstream.

The upstream `LICENSE` is preserved verbatim as [LICENSE](LICENSE). Every Skill remains under its
original licence. Per-Skill provenance — the upstream repository, the exact commit, the upstream path,
and whether the content was modified — is recorded in [PROVENANCE.json](PROVENANCE.json).

**We did not write these Skills.** What we added is the audit evidence, and the fixes that evidence
justified.

## Status

| | count |
| --- | ---: |
| Skills in this repository | **143** |
| substantively modified by us, each with a fix log | 122 |
| unmodified apart from a declared `license: MIT` | 21 |
| **fix pass still needed** (first audit only, whatever the score) | **21** |
| **re-audit still needed** (changed after their latest audit) | **21** |
| audit coverage | **100%** |

Each Skill carries a score and grade from a skill-auditor run that executed the Skill's own code against
real or synthetic data and recorded what it produced — not a review of the prose. A Skill flagged
`fix_pass: needed` still has open findings from that audit; a high score does not exempt it. A Skill
flagged `reaudit: needed` was changed after its latest audit, so its score describes earlier bytes.

Audited Skills that failed (open P0 or not deployable) are **not here**: they are listed in
`REMAINING.json` under `excluded`, with their scores, until a fix passes re-audit.

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
