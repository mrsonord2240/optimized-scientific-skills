# Internal-baseline inference
~~~sh
Rscript scripts/copykat_run.R counts.tsv tumor1 copykat_result.rds
~~~
Counts are raw genes by cells with gene symbols. Source settings use ngene.chr=5, win.size=25, KS.cut=0.1, euclidean distance, human genome hg20 and four cores. Internal baseline can fail when normal cells are absent or tumors are CNV-quiet. For SCEVAN:
~~~sh
Rscript scripts/scevan_run.R counts.tsv tumor1 scevan_result.rds normals.txt
~~~
The optional normals file has one confident normal barcode per line; omit it for internal selection. SCEVAN writes additional output under its working directory. Errors are surfaced; a console count alone is not a recovered per-cell result. Match species, gene identifiers and breadth; verify APIs.
**Done when:** returned predictions, filtered cells and baseline assumptions are checked and independently corroborated.
