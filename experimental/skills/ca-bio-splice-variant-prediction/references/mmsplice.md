## MMSplice for Calibrated ΔPSI

**Goal:** Predict quantitative ΔPSI (not just probability of disruption) for cassette exons.

**Approach:** Score variant impact on each splicing region (5'ss, 3'ss, exon, intron-3'/5') and combine.

```python
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save

dl = SplicingVCFDataloader(
    gtf='gencode.v45.basic.annotation.gtf',
    fasta_file='GRCh38.primary_assembly.genome.fa',
    vcf_file='input.vcf'
)

model = MMSplice()
predict_save(model, dl, 'mmsplice_predictions.csv', pathogenicity=True)
```

MMSplice (Cheng 2019 *Genome Biol*) reports Δlogit_psi per variant. Useful when calibrated effect sizes matter (research) more than probability of disruption (clinical screening). The CSV has **one row per variant x exon x transcript**, ID `chrom:pos:ref>alt`, and no row for a variant that is not near an annotated exon (the deep-intronic GLA c.639+919G>A got none, silently): reduce with `parse_mmsplice_csv()` (largest |delta_logit_psi| per variant) and report variants missing from the output. A truncated or mismatched FASTA yields empty sequences and assertion errors (`Input sequence ... intron length cannot be longer than the input sequence`). Companion **MTSplice** (Cheng 2021 *Genome Biol*) adds tissue-specific Δψ predictions.
