# Percolator picked proteins
~~~bash
percolator -f target_decoy.fasta -P DECOY_ -z trypsin -l prot.target.tsv -L prot.decoy.tsv -r pep.target.tsv -B pep.decoy.tsv -m psm.target.tsv -M psm.decoy.tsv -S 1 search.pin
~~~
Match FASTA, decoy prefix and enzyme to the search; verify nonempty result files and relevant group-member options. The source also provides a conditional Philosopher sequence in [external methods](../references/percolator-philosopher-cli.md); check each output rather than exit status alone.
Done when output contents, group representation and FDR context are reviewed.
