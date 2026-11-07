# VCF normalization versus HGVS

No script. VCF left-alignment moves an ambiguous indel to the most 5' position on the forward genomic strand. HGVS requires the 3'-rule: the most 3' position in the transcript's reading direction. For a plus-strand gene those are opposite ends of the repeat; for a minus-strand gene they can coincide by accident of strand. So a correctly normalized VCF and a correct HGVS string can legitimately name different positions.

- Never derive HGVS from POS. Left-align the VCF for matching and merging, and let the annotation engine produce HGVS (VEP applies the transcript 3'-shift; verify SnpEff and ANNOVAR for your version).
- To match a patient variant to a ClinVar or literature HGVS string, normalize both to the same representation first; two strings that look different can be one variant.

Done when HGVS strings come from the annotation engine, and any comparison used normalized forms on both sides.
