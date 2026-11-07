# Count reads
~~~bash
mageck count --list-seq library.csv --sample-label Plasmid,Day0,Veh_r1,Veh_r2,Drug_r1,Drug_r2 --fastq Plasmid.fq.gz Day0.fq.gz Veh_r1.fq.gz Veh_r2.fq.gz Drug_r1.fq.gz Drug_r2.fq.gz --norm-method median --output-prefix screen --trim-5 AUTO
~~~
Match each label to a FASTQ and provide library columns in guide ID, sequence, gene order with a header. Inspect countsummary, mapping, zero counts and distribution; see canonical bio-crispr-screens-screen-qc. Under heavy selection consider --norm-method control with --control-sgrna ntcs.txt. The source's version is MAGeCK 0.5.9+; confirm local help.
Done when counts, sample mapping, normalization and QC are reviewed.
