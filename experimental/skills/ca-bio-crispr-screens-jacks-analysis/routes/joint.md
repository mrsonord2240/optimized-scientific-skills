# Fit joint effects
~~~bash
python scripts/run_jacks_joint.py counts.txt replicatemap.txt guidemap.txt --outprefix jacks_out
~~~
Counts are guide rows by sample; replicate map has Replicate, Sample, Control and guide map has sgRNA, Gene headers. Match sample and guide identifiers and library and chemistry. This wrapper leaves the hierarchical gene-effect prior off and defaults to zero pseudo-genes. For p-values supply --ctrl-genes NEGv1.txt --n-pseudo 2000 --seed 1 and set PYTHONHASHSEED=1 before Python starts. Check convergence, missing guides and installed JACKS 0.2 API. Outputs include gene effects, posterior std and guide efficacy. See [model and interpretation](../references/model.md).
Done when input mapping, convergence, reproducibility and effect/std interpretation are reviewed.
