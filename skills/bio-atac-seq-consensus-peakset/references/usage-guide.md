# Consensus peakset request examples

Use these concise requests to specify the desired output. The core Skill describes the shared workflow; [`method-reference.md`](method-reference.md) contains method details and limitations.

- "Build a Corces 2018 iterative-overlap consensus peakset (501 bp fixed-width) from per-replicate narrowPeak files."
- "Use DiffBind `dba.count(summits=250)` for fixed-width counting."
- "Build a consensus per condition, then union across conditions to retain condition-specific peaks."
- "Build an IDR-filtered union from true replicate pairs at threshold 0.05."
- "Filter the consensus against the assembly-matched blacklist and convert BED to SAF for featureCounts."
- "Compare peaks across genome builds; first confirm whether the published peaks can be lifted over and document unmapped regions."

Provide the input peak files, sample-to-condition mapping, genome assembly, chromosome sizes, intended downstream analysis, and desired output format. Ask for clarification when those choices are missing and change the strategy.
