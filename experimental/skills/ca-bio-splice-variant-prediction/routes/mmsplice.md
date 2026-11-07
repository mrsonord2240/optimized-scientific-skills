# Exon-focused quantitative prediction

```python
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_save
dl = SplicingVCFDataloader(gtf='annotation.gtf', fasta_file='genome.fa', vcf_file='input.vcf')
predict_save(MMSplice(), dl, 'mmsplice_predictions.csv', pathogenicity=True)
```

Use a matching GTF and FASTA; compare output IDs with input because variants far from an annotated exon can be skipped. See `references/mmsplice.md` for output interpretation.

Done when requested predictions and unscored variants are identified.
