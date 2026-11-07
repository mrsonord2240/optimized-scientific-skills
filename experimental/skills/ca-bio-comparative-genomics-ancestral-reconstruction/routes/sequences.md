# Infer sequence states
```bash
bash scripts/iqtree_ancestral.sh
```
**Done when:** the output .state table and tree are available and node/site probabilities, rooting and model assumptions have been examined. The source script uses fixed filenames `alignment.fasta`, `outgroup_taxon`, model selection and prefix; edit them for the task before running and verify the installed IQ-TREE interface. Compare plausible alignments, roots and models; a best state is an estimate, not an observed ancestor. IQ-TREE posterior parsing is in `scripts/iqtree_state.py`.
For the source GRASP example use `bash scripts/grasp_asr.sh alignment.fasta species.nwk output-dir 4` after checking its interface and expected outputs; it validates output artifacts. PAML codeml helpers are functions requiring task-specific control-file setup and invocation; see [method workflows](../references/method-workflows.md). The PAML parser reports the best state and probability, not the full posterior distribution.
