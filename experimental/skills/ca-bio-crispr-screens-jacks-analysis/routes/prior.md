# Reuse an efficacy prior
~~~bash
python run_JACKS.py counts.txt replicatemap.txt guidemap.txt --rep_hdr Replicate --sample_hdr Sample --ctrl_sample_hdr Control --sgrna_hdr sgRNA --gene_hdr Gene --reffile reference_grna_JACKS_results.txt --outprefix jacks_out
~~~
Run the installed JACKS CLI from its JACKS/jacks directory, using valid paths for files. The reference grna table must have sgrna, X1 and X2 and match library sequences and chemistry, not merely reused IDs. The bundled joint wrapper does not expose --reffile. Check the installed version and reference coverage; see [prior and diagnostics](../references/efficacy-prior-and-diagnostics.md).
Done when prior identity, settings, convergence and results are documented and reviewed.
