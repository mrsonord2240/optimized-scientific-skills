# Position-only overlap

## Finding Shared Positions in Python (cyvcf2)

For a quick position-level overlap count without writing files, cyvcf2 is convenient (for true set operations on normalized inputs, prefer `bcftools isec`):

```python
from cyvcf2 import VCF

vcf1_positions = {(v.CHROM, v.POS) for v in VCF('sample1.vcf.gz')}

shared = unique = 0
for v in VCF('sample2.vcf.gz'):
    if (v.CHROM, v.POS) in vcf1_positions:
        shared += 1
    else:
        unique += 1

print(f'Shared positions: {shared}')
print(f'Unique to sample2: {unique}')
```

Position-only overlap ignores REF/ALT, so normalize both inputs first if allele-level concordance is the question.
