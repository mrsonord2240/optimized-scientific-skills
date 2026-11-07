## NMD Prediction (The 50-nt Rule)

A transcript is predicted NMD-sensitive if its premature termination codon (PTC) lies **>50-55 nt upstream of the last exon-exon junction** (Maquat 2004 *Nat Rev Mol Cell Biol*; Lykke-Andersen & Jensen 2015 *Nat Rev Mol Cell Biol*).

**Mechanism:** Spliceosome deposits the Exon Junction Complex (EJC) ~20-24 nt upstream of every exon-exon junction. During the pioneer round of translation, ribosome reading through removes EJCs upstream of the stop codon. If a stop codon precedes the last EJC by >50 nt, the EJC remains, recruits UPF1 -> SMG1 phosphorylation -> SMG6/SMG7 -> mRNA decay.

**Caveats and exceptions:**
- **Last-exon PTCs escape NMD** — can be dominant-negative or gain-of-function (e.g. MYH7 truncating variants).
- **3'UTR length matters**: very long 3' UTRs (>1 kb past stop) trigger NMD via UPF1 binding even without EJCs (faux-3'UTR rule).
- **Tissue-specific NMD**: SMG6 vs SMG5/7 ratios vary; UPF1 stress conditions modulate.
- **PTC distance must be measured on the spliced transcript**, not the genomic distance.
- **Some predicted-NMD transcripts escape decay in RNA-seq** (Lindeboom 2016 *Nat Genet*; the escape fraction is tissue- and rule-dependent and was not verified here). Treat NMD prediction as probabilistic, not certain.

`analyzeSwitchConsequences` with `'NMD_status'` evaluates this from the ORF (annotated if the GTF has CDS, else predicted; see the ORF order above) and the transcript model. Checked: in the shipped demo the PTC flag marks exactly the 8 planted poison isoforms, and on real chrX the annotated-ORF PTC flag agrees with Ensembl's NMD biotype (sensitivity 0.96, above).

## AS-NMD as a Regulatory Layer

A large class of conserved alternative splicing events is **deliberately PTC-introducing** to titrate functional protein levels:

- **All major SR proteins** (SRSF1-12) autoregulate via poison exons (Lareau 2007 *Nature*; Ni 2007 *Genes Dev*)
- **All major hnRNPs** likewise
- **Ribosomal protein genes** use AS-NMD autoregulation (e.g. rpL3, rpL12; Cuccurese 2005 *NAR*)
- **SCN1A** poison exon -> STK-001 (Stoke Therapeutics) ASO for Dravet syndrome, in clinical development (Han 2020 *Sci Transl Med*; trial phase not verified here)

**Functional implication:** an *increase* in PSI of a poison exon *decreases* functional protein. Sign-of-effect in DTU output is opposite from intuition for these genes. Always check whether the alternative form is PTC-bearing before interpreting direction.

**Disease examples:**
- TDP-43 cryptic exon in UNC13A introduces a PTC -> NMD of a disease-relevant transcript (Brown 2022 *Nature*). The STMN2 cryptic exon works differently: premature polyadenylation gives a truncated transcript, not PTC-NMD (Melamed 2019 *Nat Neurosci*; Klim 2019 *Nat Neurosci*). These appear in post-mortem ALS brain, not blood: check the tissue.
- Last-exon truncating variants in TTN (and MYH7): escape NMD -> stable poison/dominant-negative protein
