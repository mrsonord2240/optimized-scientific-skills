# Review substitution saturation
~~~bash
Rscript scripts/saturation.R alignment.fasta saturation_pairs.tsv
~~~
The source ape proxy uses raw p-distance, transitions and transversions divided by alignment length, and transition slope versus TN93. Review missing/nonfinite distances and missing-data handling; this is not Xia's formal test. A finite corrected distance or strong bootstrap does not establish absence of saturation. See [saturation guidance](../references/methods.md).
Done when the diagnostic, model assumptions and limitations are documented.
