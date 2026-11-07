# Identity-distance neighbor joining
~~~bash
python scripts/nj_identity.py alignment.fasta tree.nwk
~~~
The source procedure removes every gap-containing alignment column, computes Biopython identity distance and builds NJ. Identity distance is uncorrected and is p-distance only under suitable gap/ambiguity handling; it is not JC/K80. Inspect retained sites and sequence quality.
Done when matrix, gap policy and tree assumptions are recorded.
