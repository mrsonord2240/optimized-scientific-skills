# cyvcf2 recipes

## cyvcf2 Python Alternative

**Goal:** Read, query, and write VCF files programmatically in Python.

**Approach:** Use cyvcf2's `VCF` reader to iterate variants with attribute access to fields, and `Writer` to emit filtered output.

**"Parse this VCF in Python"** -> Open with cyvcf2 and iterate variant records.

### Open, Iterate, and Access Fields
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
for variant in vcf:
    # ALT is a list; QUAL is site-level and may be None
    print(variant.CHROM, variant.POS, variant.REF, variant.ALT)
    print(variant.ID, variant.QUAL, variant.FILTER, variant.var_type)
    dp = variant.INFO.get('DP')   # INFO field, None if absent
    af = variant.INFO.get('AF')
    break
vcf.close()
```

### Access Genotypes and Per-Sample Fields
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
samples = vcf.samples
for variant in vcf:
    # gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN(missing), 3=HOM_ALT
    gts = variant.gt_types
    depths = variant.format('DP')   # numpy array, one row per sample
    gqs = variant.format('GQ')      # per-sample genotype quality
    ad = variant.format('AD')       # per-allele depth, Number=R
    print(dict(zip(samples, gts)))
    break
vcf.close()
```

Note: cyvcf2 codes missing genotypes as `gt_types == 2` (UNKNOWN) -- treat that as no-call, never as HOM_REF.

### Fetch Region and Read the Header
```python
from cyvcf2 import VCF

vcf = VCF('input.vcf.gz')
print(vcf.samples, vcf.seqnames)      # sample names, contig names
for info in vcf.header_iter():
    if info['HeaderType'] == 'INFO':
        print(info['ID'], info['Description'])
for variant in vcf('chr1:1000000-2000000'):   # requires an index
    print(variant.CHROM, variant.POS)
```

### Write Filtered VCF
```python
from cyvcf2 import VCF, Writer

vcf = VCF('input.vcf.gz')
writer = Writer('output.vcf', vcf)   # inherit the input header
for variant in vcf:
    if variant.QUAL is not None and variant.QUAL > 30:   # QUAL is site-level
        writer.write_record(variant)
writer.close()
vcf.close()
```


## References

- Danecek P, Auton A, Abecasis G, et al. The variant call format and VCFtools. *Bioinformatics.* 2011;27(15):2156-2158. doi:10.1093/bioinformatics/btr330 (VCF format definition)
- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008 (bcftools view/query/norm reference)
- The Variant Call Format Specification (VCFv4.3/4.4). GA4GH / samtools hts-specs. https://samtools.github.io/hts-specs/ (symbolic alleles, END, `*` overlapping-deletion allele, Number=A/R/G, PL/GL/GQ, gVCF `<NON_REF>`)
