# Base and prime editing
~~~bash
CRISPResso --fastq_r1 cbe_sample.fastq.gz --amplicon_seq "$AMPLICON" --guide_seq "$GUIDE" --base_editor_output --conversion_nuc_from C --conversion_nuc_to T --quantification_window_size 10 --quantification_window_center -10 --output_folder cbe_results --name cbe_sample
~~~
Set AMPLICON and GUIDE to the input sequences. For ABE use A to G; separately report target, bystander and indel rates. For prime editing, supply spacer, extension and scaffold sequences with the flags and interpretation in [prime editor](../references/prime-editor.md); see [base editor](../references/base-editor.md). Window center differs between source examples and its cited guidance; verify positions and installed-version semantics.
Done when intended, bystander or scaffold outcomes and the window are reviewed.
