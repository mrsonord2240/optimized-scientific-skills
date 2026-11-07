# VCF field and genotype semantics

## Format Overview

| Format | Description | Use Case |
|--------|-------------|----------|
| VCF | Text format, human-readable | Debugging, small files |
| VCF.gz | Compressed VCF (bgzip) | Standard distribution |
| BCF | Binary VCF | Fast processing, large files |

## VCF Format Structure

```
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE1
chr1    1000    rs123   A       G       30      PASS    DP=50   GT:DP   0/1:25
```

### Header Lines (##)
- `##fileformat` - VCF version
- `##INFO` / `##FORMAT` - INFO / FORMAT field definitions (ID, Number, Type)
- `##FILTER` - Filter definitions
- `##contig` - Reference contigs (required for region indexing and contig order)
- `##reference` - Reference genome

### The Header Contract

For portable VCFs, declare every INFO/FORMAT tag used in the body in a `##INFO`/`##FORMAT` line giving its ID, Number, and Type. VCF permits absent meta-information, but BCF requires its dictionaries, and parsers (bcftools, cyvcf2, pysam) use these declarations to type and reshape fields. An out-of-sync header -- especially a wrong Number/Type such as `Number=1` over vector data -- can make a tool drop, mistype, or mis-subset values without a hard error. `##contig` declarations are recommended for VCF and required for BCF; include them for indexed/convertible files and known contig order. After any hand-edit or annotation that adds a field, update the header to match (`bcftools +fill-tags` and `bcftools annotate` manage this automatically).

### Data Columns

| Column | Description |
|--------|-------------|
| CHROM | Chromosome |
| POS | 1-based position of the first base in REF (contrast BED's 0-based half-open) |
| ID | Variant identifier (e.g., rs number) or `.` if novel |
| REF | Reference allele (matches the reference exactly) |
| ALT | Alternate allele(s), comma-separated. `*` = allele missing due to an overlapping deletion at this site |
| QUAL | Phred-scaled quality of the ALT assertion, `-10*log10 P(no variant)`; higher = more confident a variant exists (site-level, NOT per-sample) |
| FILTER | PASS or semicolon-separated filter names. `.` means filters were not applied |
| INFO | Semicolon-separated key=value pairs (site-level annotations) |
| FORMAT | Colon-separated format keys defining per-sample field order |
| SAMPLE | Colon-separated values matching FORMAT order |

## Critical Field Interpretation

What each field actually measures -- and what it does not -- drives every filtering and interpretation decision. QUAL, GQ, and PL answer three DIFFERENT questions.

### QUAL vs GQ vs PL/GL: three different confidences

| Field | Level | Scale | Question answered |
|-------|-------|-------|-------------------|
| QUAL (col 6) | Site | Phred: `-10*log10 P(no variant)` | "Is there ANY variant at this site?" |
| GQ (FORMAT) | Genotype | Phred, capped at 99 | "Is THIS sample's assigned genotype correct?" |
| PL (FORMAT) | Genotype | Phred, rebased to min=0 | Relative likelihood of every possible genotype |
| GL (FORMAT) | Genotype | log10 likelihoods, usually `<=0` | Same likelihood information as PL before relative phred scaling and integer rounding |

QUAL is computed once across all samples and SCALES with total depth, so a high-coverage artifact can carry a large QUAL -- hence QD (QUAL normalized by depth) is preferred for filtering. GQ is per-sample and does not scale with cohort size. They are NOT interchangeable: QUAL can be high while an individual genotype is uncertain (low GQ), and a sample can have a confident genotype (high GQ) at a site with only moderate QUAL. Filter site-level junk on QUAL/QD; no-call untrustworthy genotypes on GQ.

### PL/GL, and how GQ is derived

PL holds phred-scaled genotype likelihoods, rebased so the best-supported genotype has 0 and every other value is its phred penalty relative to it. For a biallelic diploid site PL is ordered `[PL(0/0), PL(0/1), PL(1/1)]`; the minimum-PL index normally agrees with GT, but GT remains the caller's emitted genotype and is authoritative when they disagree. GL stores log10 likelihoods (larger is better). Given comparable GL values, PL is approximately `round(-10 * (GL - max(GL)))`; rounding means it is not simply `-10*GL`.

Many callers derive GQ from the difference between the two smallest PL values, i.e. phred separation between the best and next-best genotype; GQ=0 means a tie after rounding and it is commonly capped at 99. Treat this as a useful validation check for the caller at hand, not a universal reconstruction rule: callers may compute or cap GQ differently, and missing/partial likelihood vectors cannot establish GQ.

For a site with n alleles, diploid genotype `j/k` (j<=k) sits at PL index `k*(k+1)/2 + j` (this is the `Number=G` ordering). Getting this index formula wrong is the classic bug when re-parsing PL after a multiallelic split -- the vector must be re-subset by the formula, never sliced positionally.

### AD vs DP, and allele balance

AD (FORMAT, `Number=R`) is per-allele read depth `[ref_depth, alt1_depth, ...]`, REF first. DP is total depth. `sum(AD)` is often LESS than DP -- expected, not an error:
- DP counts all reads spanning the position, including uninformative reads (low base quality, ambiguous alignment, filtered reads).
- AD counts only reads that confidently support a specific allele.
- INFO/DP (site-level, summed across samples) differs from FORMAT/DP (per-sample).

Allele balance for a het is DERIVED from AD (GATK does not emit it directly): `AB = alt_AD / (ref_AD + alt_AD)`. A true het sits near 0.5; hets far from 0.5 (e.g. `<0.2` or `>0.8`) suggest a mapping artifact, CNV, or contamination.

### INFO/FORMAT Number semantics (A / R / G / .)

Every `##INFO`/`##FORMAT` header declares a `Number` telling a parser how many values a field holds AND how to re-subset it when a multiallelic record is split:

| Number | One value per | Examples | On multiallelic split |
|--------|---------------|----------|-----------------------|
| `A` | ALT allele | AF, AC | take the k-th element for the k-th ALT |
| `R` | allele incl. REF | AD | REF value first, then per-ALT (off-by-one vs A) |
| `G` | genotype | PL, GL | re-subset via the `k*(k+1)/2+j` index formula |
| `.` | variable/unknown | -- | parser CANNOT auto-subset; carried whole onto every record |
| `0` | flag (presence only) | -- | -- |

Load-bearing for correctness: `bcftools norm -m-` uses these codes to reapportion fields on split. A field mis-declared `Number=.` when it is really `A` keeps its full multiallelic vector on every split record, so downstream tools read the WRONG allele's value with no error. See variant-calling/variant-normalization for split/join reapportionment.

### Key INFO Annotations for Filtering

| Annotation | Meaning | What It Detects |
|-----------|---------|-----------------|
| QD | QUAL / allele depth | Low values suggest variant quality not supported by reads |
| FS | Fisher strand bias (phred-scaled) | Variant reads predominantly on one strand (artifact) |
| SOR | Strand odds ratio | Same as FS but handles high-depth sites better |
| MQ | Root mean square mapping quality | Low values indicate reads map ambiguously (paralogous regions) |
| MQRankSum | MQ difference: ref vs alt reads | Very negative = alt reads map much worse than ref (suspicious) |
| ReadPosRankSum | Read position: ref vs alt reads | Very negative = variant only at read ends (misalignment artifact) |

## Genotype Encoding

| Genotype | Meaning |
|----------|---------|
| `0/0` | Homozygous reference (confidently called ref) |
| `0/1` | Heterozygous |
| `1/1` | Homozygous alternate |
| `1/2` | Heterozygous for two different ALT alleles (compound het at multiallelic site) |
| `./.` | Missing genotype (no confident call) |
| `0\|1` | Phased heterozygous (allele before `\|` is on haplotype 1) |

### Phased vs Unphased

- `/` separates **unphased** alleles -- the two chromosomal copies are known, but which came from which parent is not
- `|` separates **phased** alleles -- haplotype assignment is known (read-backed phasing, trio analysis, or long-read sequencing)
- Phasing matters for compound heterozygosity: two variants in a gene are pathogenic together only if on different haplotypes (in *trans*), not the same haplotype (in *cis*)

### Phase Sets (PS)

A `|` is only meaningful WITHIN a phase set. The FORMAT/PS tag (an integer, usually the POS of the block's first variant) groups variants phased relative to EACH OTHER; `0|1` in two different PS blocks are not guaranteed to lie on the same physical haplotype. Read-backed phasers (WhatsHap) and trio phasing emit PS. A `|` with no consistent PS across records carries no global phase -- a subtle trap when merging phased VCFs.

### Ploidy and Missing vs Hom-Ref

Ploidy is read from the NUMBER of alleles in GT: `0/1` diploid, `0` haploid (chrY, chrM, male chrX outside the PAR), `0/1/1` triploid. Per-region ploidy (PAR, chrX in males, mito) must match the sample karyotype.

`.` (missing) is NOT reference. `./.` = no-call (genotype could not be determined, usually low depth); `0/0` = confidently called homozygous reference. Treating `./.` as `0/0` inflates the reference-allele count and biases allele frequencies, missingness, and burden tests. This is load-bearing: never impute `./.` as reference. In a gVCF, the ABSENCE of a record also does not mean reference -- see the reference-confidence model below.

### Multiallelic Genotypes

At multiallelic sites (e.g. ALT = G,T), allele indices reference the comma-separated ALT list: 0=REF, 1=first ALT, 2=second ALT. `1/2` means one copy of each ALT. Splitting multiallelics into biallelic records with `bcftools norm -m-` converts `1/2` into two `0/1` records, losing compound-heterozygosity information -- see variant-calling/variant-normalization for caveats.

## Symbolic Alleles, END, and Spanning Deletions

Not every ALT spells out a sequence. Symbolic alleles are angle-bracketed placeholders for events whose sequence is not given inline:

| ALT | Meaning |
|-----|---------|
| `<DEL>` `<DUP>` `<INS>` `<INV>` `<CNV>` | Structural-variant classes (sequence not spelled out) |
| `<NON_REF>` | gVCF: "any allele not yet observed" (reference-confidence model) |
| `<*>` | Same role as `<NON_REF>` in some callers' gVCF/mpileup output |
| `*` (bare) | Spanning deletion: allele MISSING because an upstream deletion on ANOTHER line overlaps this position |

Two parsing traps:
- `INFO/END` gives the end coordinate of a symbolic/large event. A tool that infers a record's span from `len(REF)` is WRONG for symbolic alleles -- it must read END. gVCF reference blocks also use END to mark the last position of the band.
- The bare `*` ALT is interpretable only relative to the overlapping deletion on another record; it is not a real alternate allele here. Splitting/subsetting can strand a `*` from the deletion it references (see variant-calling/variant-normalization).

### gVCF and the `<NON_REF>` Reference-Confidence Model

A gVCF (GATK HaplotypeCaller `-ERC GVCF`) is fundamentally different from a filtered callset: it emits a record for EVERY position or block, not just variant sites. Non-variant stretches are compressed into END-delimited blocks (bands) grouped by GQ, so a gVCF is not one line per base.

- Every record carries a symbolic `<NON_REF>` ALT with PL/AD computed against "any unseen allele." This lets joint genotyping evaluate a site in THIS sample even when the variant was only discovered in ANOTHER cohort sample -- the `<NON_REF>` likelihood supplies the evidence.
- Its purpose is to distinguish, at every site, confident homozygous reference from no-data/no-call -- solving the missing-vs-reference problem when squaring off a cohort matrix.
- A gVCF is NOT ready for analysis; it is an intermediate. It must be joint-genotyped (`GenomicsDBImport`/`CombineGVCFs` -> `GenotypeGVCFs`) to yield a normal VCF. Do NOT filter, annotate, or count variants on a raw gVCF, and never build a multi-sample callset by `bcftools merge`-ing single-sample project VCFs when gVCF joint-genotyping is available -- merging fabricates hom-ref genotypes. See variant-calling/joint-calling.
- To list candidate variant sites (not reference bands) before joint genotyping: `bcftools view -i 'N_ALT>1' sample.g.vcf` -- a record has more than one ALT only when a real allele sits alongside `<NON_REF>`. A string test like `ALT="<NON_REF>"` is NOT a substitute: it matches on presence in the ALT list, so it matches every record, reference bands included, not just candidate sites.


## References

- Danecek P, Auton A, Abecasis G, et al. The variant call format and VCFtools. *Bioinformatics.* 2011;27(15):2156-2158. doi:10.1093/bioinformatics/btr330 (VCF format definition)
- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008 (bcftools view/query/norm reference)
- The Variant Call Format Specification (VCFv4.3/4.4). GA4GH / samtools hts-specs. https://samtools.github.io/hts-specs/ (symbolic alleles, END, `*` overlapping-deletion allele, Number=A/R/G, PL/GL/GQ, gVCF `<NON_REF>`)
