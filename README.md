# Optimized Scientific Skills

Agent Skills for computational science that have been **audited against executed evidence**. Where the
audit demonstrated a defect, the Skill is fixed and then re-audited by an agent that did not write the fix.

## Project status

The audit/evaluation campaign was closed out on **2026-11-07**, preserving the
195 audited scientific Skills, six process Skills, and evaluation evidence in
the [pipeline archive](pipeline/README.md).

Sam then authorized publishing the **60 partially complete normalized
candidates** under [experimental/](experimental/README.md), and resumed
normalization of the other 171 Skills already on this shelf. These drafts have
their own status notes; their changed bytes do not inherit the published
versions' audit grades. No new scientific audit or Marketplace acceptance is
claimed. [PROVENANCE.json](PROVENANCE.json) describes the 195 audited versions.

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

These counts describe the 195 audited versions under `skills/`, excluding the experimental drafts:

| | count |
| --- | ---: |
| Audited Skills under `skills/` | **195** |
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
experimental/             partially complete ca- revisions, individual status notes and source hashes
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

## Experimental revisions

The [60 experimental candidates](experimental/README.md) are now included with
their recorded limitations and individual status notes. New normalization-only
work is published in the same section as it is reviewed. The evaluation's draft
arms also remain experimental. The [evidence index](pipeline/EVIDENCE.md)
records the original source revisions and local preservation archives.
