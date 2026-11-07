# Optimized Scientific Skills

Agent Skills for computational science that have been **audited against executed evidence**. Where the
audit demonstrated a defect, the Skill is fixed and then re-audited by an agent that did not write the fix.

## Project status

This work was closed out on **2026-11-07**. The 195 scientific Skills and their
evidence remain available. The preserved [pipeline](pipeline/README.md) contains
the six process Skills, the [evaluation method and results](pipeline/EVALUATION.md),
and an [evidence and recovery index](pipeline/EVIDENCE.md). Further audits,
experimental recuts, and Open Science Marketplace follow-up are no longer an active queue.

The closeout preserves prior results; it is not a new scientific audit or
Marketplace acceptance. Each Skill's recorded status is in
[PROVENANCE.json](PROVENANCE.json). [REMAINING.md](REMAINING.md) retains the
historical unprocessed/excluded inventory, not a commitment to resume it.

## Attribution

These Skills originate in **[GPTomics/bioSkills](https://github.com/GPTomics/bioSkills)**, MIT
licensed, by GPTomics. That repository was archived on 2026-08-15 and accepts no issues or pull
requests, which is why refinement happens here rather than upstream.

The upstream `LICENSE` is preserved verbatim as [LICENSE](LICENSE). Every Skill remains under its
original licence. Per-Skill provenance — the upstream repository, the exact commit, the upstream path,
and whether the content was modified — is recorded in [PROVENANCE.json](PROVENANCE.json).

The upstream authors retain credit for the scientific Skills. Our contribution
is refinement, audit evidence, and process/evaluation tooling. The separately
preserved pipeline has its own [license and attribution](pipeline/NOTICE).

## Status

These counts describe the committed provenance snapshot, not fresh execution:

| | count |
| --- | ---: |
| Skills in this repository | **195** |
| marked modified, each with a fix-log reference | 195 |
| unmodified apart from declared `license: MIT` and `author: GPTomics` lines | 0 |
| **fix pass still needed** (first audit only, whatever the score) | **0** |
| **re-audit still needed** (changed after their latest audit) | **0** |
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
pipeline/                 preserved process Skills, evaluation, and evidence/recovery index
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

## Preserved experiments

The records repository retains 24 `ca-` experimental recuts with their individual
states. They are not installed into this repository's scientific Skill collection.
The evaluation's draft arms also remain experimental. The
[evidence index](pipeline/EVIDENCE.md) pins their source revision and records the
local preservation archives. This closeout does not change scientific Skill bytes
or their provenance.
