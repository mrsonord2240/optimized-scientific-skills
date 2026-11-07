# Transcript set, `--pick`, and HGVS

No script. Decide, then record the choice with the result.

| Set | Use |
|-----|-----|
| MANE Select | Default for clinical reporting on GRCh38: one transcript per gene, identical `NM_` and `ENST`, stable c./p. |
| MANE Plus Clinical | Add alongside MANE Select where assigned; it covers known pathogenic variants MANE Select misses |
| RefSeq (`NM_`) | Legacy clinical pins; most ClinVar and literature c. notation |
| Ensembl/GENCODE (`ENST`) | Research; not 1:1 with RefSeq, c. positions differ. Never compare a RefSeq-annotated variant with an Ensembl one |
| Worst consequence across all | Discovery only: inflates severity and manufactures false PVS1 candidates |

- The Ensembl canonical transcript is often not the MANE Select one, so moving a pipeline to MANE changes some reported c. coordinates. That is expected. MANE is GRCh38 only; GRCh37 pipelines lift over or keep a per-gene transcript pin.
- `--pick` collapses to one block per variant by an ordered heuristic. Since VEP 110 the default `--pick_order` leads with `mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq`; older releases fell through to canonical, biotype, rank, then length and accession. Late tiebreakers let transcript length or `ENST` order choose a patient's consequence, and different variants in one gene can land on different transcripts. If used, constrain it: `--mane_select --pick --pick_order mane_select,mane_plus_clinical,canonical,biotype,rank`. Or pin MANE (with Plus Clinical), or the lab's validated per-gene list.
- VCF left-alignment and HGVS conflict: VCF puts an indel in a repeat at the most 5' genomic position; HGVS puts it at the most 3' position on the transcript. A correct POS and a correct `c.` can then point to different repeat units. Left-align before annotating, let the engine apply the transcript 3'-shift (VEP does; verify SnpEff and ANNOVAR per version), and never hand-derive HGVS from POS. To match against ClinVar or literature, normalize both sides; `dup` and `ins` describe one event but do not string-match.
- Adjacent SNVs in one codon split by the caller can each look benign and jointly change the amino acid; use phase-aware annotation.

Done when the report names the transcript set, the selection rule, and the engine and version.
