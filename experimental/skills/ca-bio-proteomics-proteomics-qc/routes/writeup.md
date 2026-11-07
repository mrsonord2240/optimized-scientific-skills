# QC write-up

Fill `references/qc-report-template.md`: metric table (what was and was not measurable), exclusion and decision log, with/without sensitivity table, design limits, result.

- One log row per flagged sample or channel: value, threshold, decision (exclude, re-run, keep-flagged, keep), who approved. Exclude loading failures, never rescale them.
- Run the downstream check with and without every excluded or keep-flagged sample and report both. If the conclusion changes, say so.
- State design limits in words: groups of one, a batch confounded with condition, fewer than 5 per group, no un-normalized column.
- Thresholds quoted in the write-up: `references/thresholds.md`.

Done when the write-up has all six sections and no unmeasurable check is reported as passing.
