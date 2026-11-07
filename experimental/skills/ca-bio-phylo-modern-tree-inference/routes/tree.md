# Infer an ML tree
~~~bash
iqtree3 -s alignment.fasta -m MFP -B 1000 -bnni -alrt 1000 -T AUTO -ntmax 8 --seed 12345 --prefix run1
~~~
Use an inspected alignment and a fresh output prefix. The source uses ModelFinder, 1000 UFBoot and SH-aLRT replicates with NNI refinement; select installed IQ-TREE executable and verify flags. Interpret joint SH-aLRT >=80 and UFBoot >=95 as source conventions, not proof of correctness. Inspect informative sites and near-zero branches. For RAxML-NG alternative see [methods](../references/methods.md).
Done when input, selected model, support, logs and model sensitivity are reviewed.
