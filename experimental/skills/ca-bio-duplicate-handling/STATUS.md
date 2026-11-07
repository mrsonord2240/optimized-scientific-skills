# Experimental status: ca-bio-duplicate-handling

Published as partially complete work on 2026-11-07. This is a context-aware
revision of `bio-duplicate-handling`. The earlier published Skill's audit score does not
certify these rewritten bytes. Do not present this draft as Production Ready
or Marketplace-accepted.

## Recorded state

Normalized only, not audited. Executed: all five scripts (samtools, pysam, umi_tools, fgbio on staged BAMs). Unexecuted: `umi_tools_dedup.sh scrna` success path, alternative-marker commands.

Publication preserves the recorded state; it does not rerun the scripts,
complete the missing routes or tests, or clear the listed defects. Historical
execution or evaluation claims above apply only to the described surfaces.
Tooling preparation and scientific audit are still required for unverified
surfaces before relying on their results.

## Source and credit

Preserved from [the records repository](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/skills/ca-bio-duplicate-handling). Original scientific
Skill: GPTomics/bioSkills, MIT, upstream path `alignment-files/duplicate-handling`
at `d91ed3d563019e649dc854c56ccd62551359488a`. The candidate is a modified
derivative; GPTomics retains original author credit. See the repository
[MIT license](../../../LICENSE) and [candidate manifest](../../CANDIDATES.json).
