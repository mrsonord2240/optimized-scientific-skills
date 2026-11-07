# Infer groups with OpenMS
~~~bash
python scripts/group_proteins_pyopenms.py peptides_with_decoys.idXML --decoy-prefix DECOY_ --fdr 0.01 --out groups.tsv
~~~
Use idXML containing posterior error probabilities and retained target and decoy hits. The script applies BasicProteinInferenceAlgorithm with indistinguishable grouping and greedy resolution, then picked FDR. Output includes all groups, including decoys; filter using is_decoy and qvalue. Verify decoy prefix, peptide evidence, version and output membership. See [basic inference](../references/pyopenms-basic-inference.md).
Done when input score orientation, group membership and q-values are reviewed.
