# Bayesian Inference with EPIFANY

### Bayesian Inference + Group FDR with EPIFANY

**Goal:** Assign calibrated protein/group posteriors and control protein-group FDR with a probability model rather than greedy parsimony.

**Approach:** EPIFANY consumes idXML whose PSMs already carry posterior error probabilities (from Percolator or IDPosteriorErrorProbability), then propagates belief over the peptide-protein graph. The TOPP tool is `Epifany`; the pyOpenMS class is `BayesianProteinInferenceAlgorithm`.

```bash
python scripts/epifany_inference.py peptides_with_pep.idXML --out groups.tsv
```

`scripts/epifany_inference.py` prints each group's accessions and posterior (higher = more likely present). The third positional argument of `inferPosteriorProbabilities` is `greedy_group_resolution` (`--greedy-group-resolution`); flipping it did not change the result on the reference idXML, so do not assume it removed subsumable proteins (see the EPIFANY row of the Tool Taxonomy in `SKILL.md`).
