# Search with a predicted library
~~~bash
bash scripts/diann_analysis.sh
~~~
The source example uses fixed input names, one missed cleavage, fixed and variable modifications, q=0.01, 8 threads, 15/15 ppm and reanalysis. Edit FASTA, output, instrument accuracy, modifications and inputs to match the dataset before use. DIA-NN 2.x source workflow generates a predicted library without raw files first, then searches inputs against it. Verify installed flags and output log; see [method notes](../references/methods.md).
Done when both stages, settings, output and q-value context are documented.
