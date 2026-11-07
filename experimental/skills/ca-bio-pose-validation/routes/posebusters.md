# PoseBusters validation

```python
from scripts.validate_poses import run_posebusters
results = run_posebusters('poses.sdf', 'receptor.pdb')
results.to_csv('posebusters.csv', index=False)
```

The helper is `scripts/validate_poses.py`. Pass ref_sdf='reference.sdf' for redock; without a reference, stereo, molecular identity and RMSD checks are absent. Verify loading and expected boolean check columns before trusting pb_valid. See `references/checks.md` and `references/interpretation.md`.

Done when per-pose checks, selected configuration and missing tests are reported.
