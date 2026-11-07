# RNA velocity
~~~sh
python scripts/velocity_run.py input_with_splicing.h5ad velocity.h5ad
~~~
Input must contain aligned spliced and unspliced layers; verify quantifier and preprocessing. Source settings use min_shared_counts=20, 2000 HVGs, 30 PCs/neighbors and deterministic velocity, with one job. This is a model-based estimate; inspect confidence, phase portraits and sensitivity to quantification and neighborhood choices. Historical source reports scVelo/numpy/pandas incompatibilities for stochastic and dynamical modes; check current installed versions rather than assuming those modes work.
**Done when:** output, fit diagnostics, quantifier, assumptions and orthogonal evidence have been reviewed.
