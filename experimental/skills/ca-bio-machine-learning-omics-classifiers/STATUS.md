# Experimental status: ca-bio-machine-learning-omics-classifiers

Published as partially complete work on 2026-11-07. This is a context-aware
revision of `bio-machine-learning-omics-classifiers`. The earlier published Skill's audit score does not
certify these rewritten bytes. Do not present this draft as Production Ready
or Marketplace-accepted.

## Recorded state

Normalized only, not audited. Executed once: `logistic_regression.py --signature` on Golub (dense path crashed, fixed, not rerun). Unexecuted: `rf_xgboost_classifier.py`, `batch_checks.py`, `calibration_check.py`; SMOTE path untestable (no imblearn).

Publication preserves the recorded state; it does not rerun the scripts,
complete the missing routes or tests, or clear the listed defects. Historical
execution or evaluation claims above apply only to the described surfaces.
Tooling preparation and scientific audit are still required for unverified
surfaces before relying on their results.

## Source and credit

Preserved from [the records repository](https://github.com/mrsonord2240/optimizing-agent-science-skills/tree/e9dcbbb67b1e3f02922e5d84b6778d757dce08d1/skills/ca-bio-machine-learning-omics-classifiers). Original scientific
Skill: GPTomics/bioSkills, MIT, upstream path `machine-learning/omics-classifiers`
at `d91ed3d563019e649dc854c56ccd62551359488a`. The candidate is a modified
derivative; GPTomics retains original author credit. See the repository
[MIT license](../../../LICENSE) and [candidate manifest](../../CANDIDATES.json).
