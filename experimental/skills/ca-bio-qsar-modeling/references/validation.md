## Scaffold-Balanced Split

**Goal:** Partition a SMILES dataset into train/val/test such that no Bemis-Murcko scaffold appears in more than one split (prevents chemotype leakage).

**Approach:** Group compounds by scaffold and assign whole scaffold groups to train, validation, or test. chemprop's `--split scaffold_balanced` implements a scaffold-based allocation; `--class-balance` is a separate training option and does not make this split outcome-stratified. The chemprop default split is random, so request scaffold-balanced explicitly when it matches the deployment question.

`scaffold_balanced` assigns each scaffold group to one of train / validation / test, reducing direct scaffold leakage. It is not universally the correct validation design: time splits, externally defined series, grouped cross-validation, and prospective tests may better represent a particular deployment setting.

Choose and document the primary split before model selection. A random split can answer an interpolation question but often shares close analogues across partitions; a scaffold split tests transfer across scaffold groups; a time or prospective split tests the historical deployment process. If several splits are reported, interpret their differences as split-specific sensitivity rather than a universal "true generalization gap."

This split choice also constrains conformal interpretation: scaffold and time partitions can make calibration and test examples non-exchangeable, so nominal conformal coverage is not guaranteed. Report empirical coverage on the held-out deployment-relevant partition beside the nominal level; consider Mondrian or group-conditional conformal variants when group coverage is required.
