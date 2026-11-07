# Preserved BioSkills pipeline

This work is closed out as of 2026-11-07. Sam reported no further response
from Open Science and decided the project had done enough. The retained
deliverables are the published Skills, the pipeline Skills, the audit evidence,
and the evaluation method and results. There is no active audit, recut, or
Marketplace submission queue.

## Contents

| Path | Purpose |
| --- | --- |
| [skills/](skills/) | Six original process Skills, their shared contracts, and the normalization example |
| [EVALUATION.md](EVALUATION.md) | What the evaluation actually tested, results, limits, and reproduction steps |
| [evaluation/](evaluation/) | Historical harness, synthetic task definitions/generators/solvers, and saved result records |
| [EVIDENCE.md](EVIDENCE.md) | Pinned audit records, experimental recuts, and local preservation archives |
| [MANIFEST.json](MANIFEST.json) | Source revisions and SHA-256 identities of the preserved source files |
| [LICENSE](LICENSE), [NOTICE](NOTICE) | Process/evaluation attribution and third-party license boundaries |

The process Skills are a historical snapshot from
[`optimizing-agent-science-skills@e9dcbbb67b1e3f02922e5d84b6778d757dce08d1`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1).
They retain their original text, including machine paths, phase delegation,
Marketplace gates, and publication instructions. They are stored outside the
agent discovery directories. Reading this archive does not restart that process
or authorize those historical remote actions.

## The method worth keeping

1. Normalize a Skill into a short entry point, explicit routes, runnable
   scripts, and progressively deeper references. Preserve attribution.
2. Prepare the tools and bounded inputs; record versions, environment identity,
   and surfaces that cannot be exercised.
3. Run a diagnostic audit against exact source bytes. Check meaningful outputs,
   not just imports, exit codes, or file existence.
4. Fix concrete findings. Refresh tooling only when a change affects it.
5. Have an independent auditor rerun corrected findings, likely regressions,
   and a representative core workflow. Reuse prior evidence only when its
   source and environment identity still apply.
6. Keep the report, scripts, assertions, finding dispositions, and exact source
   identity together. Candidate readiness, committing product bytes, and
   Marketplace acceptance are separate states.

The records repository contains the deterministic tooling used by these Skills:
[`tools/`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/tools),
[`scripts/`](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/scripts),
and the original `skill-auditor.zip`. They are preserved in the records Git
bundle described in [EVIDENCE.md](EVIDENCE.md), rather than duplicated here.

No new scientific audit, model evaluation, or promotion was performed for this
closeout. The 24 `ca-` recuts remain experimental in the records repository;
they are not additions to the 195 published Skills under the root `skills/`.

To inspect the saved evaluation without Docker, dependencies, or model calls:

```text
python pipeline/evaluation/summarize.py --check
```
