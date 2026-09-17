# T-Coffee: Consistency-Based Alignment Modes

Moved out of `SKILL.md` (bio-alignment-multiple) because the per-mode detail (M-Coffee, Expresso,
3D-Coffee, R-Coffee, Pro-Coffee, TCS reliability scoring) is needed less often than the default
T-Coffee run -- consult this file when the "Running T-Coffee" section in `SKILL.md` points you here.

| Mode | Flag | What it does |
|------|------|-------------|
| Default | (none) | T-Coffee + Lalign pairwise library |
| M-Coffee | `-mode mcoffee` | Combines libraries from MAFFT, MUSCLE, ClustalW, ProbCons, T-Coffee, etc. |
| Expresso | `-mode expresso` | PSI-BLAST searches PDB for structural templates, runs SAP structural alignment (requires internet for PSI-BLAST and PDB lookups; offline alternative: 3D-Coffee with user-supplied templates) |
| 3D-Coffee | `-mode 3dcoffee -template_file templates.txt` | User-supplied PDB templates; SAP / TM-align pairwise structural library |
| R-Coffee | `-mode rcoffee` | RNA: combines sequence alignment with consensus secondary structure (RNAplfold) |
| Pro-Coffee | `-mode procoffee` | Promoter regions: enforces position-specific TF binding-site alignment |
| Reliability | `-evaluate -output score_ascii` | TCS column reliability score for an existing alignment |

```bash
t_coffee input.fasta -output fasta_aln -outfile aligned.fasta

t_coffee input.fasta -mode mcoffee -output fasta_aln -outfile aligned.fasta

t_coffee input.fasta -mode expresso -output fasta_aln -outfile aligned.fasta

t_coffee -infile aligned.fasta -evaluate -output score_ascii > tcs_scores.ascii
```

**When to use T-Coffee**: Small datasets (<50 sequences) where maximum accuracy matters, especially when structural information (PDB templates) is available. Expresso (Armougom et al 2006 NAR) and 3D-Coffee modes (Poirot et al 2004; O'Sullivan et al 2004 JMB) substantially improve correct-column rate over sequence-only T-Coffee when structures exist; verify the latest benchmark numbers in the project documentation. Expresso requires internet access for PSI-BLAST + PDB lookups. The TCS reliability score (Chang et al 2014 MBE) flags individual columns as reliable/unreliable for downstream filtering before phylogenetics.

Armougom F, Moretti S, Poirot O, Audic S, Dumas P, Schaeli B, Keduas V, Notredame C. 2006. Expresso: automatic incorporation of structural information in multiple sequence alignments using 3D-Coffee. NAR 34:W604-W608.
