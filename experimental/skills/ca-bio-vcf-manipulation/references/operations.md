# Combination and subsetting details

## The governing principle: normalize BEFORE combining

Every combine operation here -- merge, concat `-d` dedup, isec, and downstream annotate -- keys on the `(CHROM, POS, REF, ALT)` tuple, and `bcftools isec` defaults to `-c none` (an ALT must match exactly to count as the same variant). An indel that is not left-aligned + parsimonious, an un-split multiallelic, or an un-decomposed MNP is a **structurally valid** VCF line carrying a *different* tuple for the *same* biological event. The combine then silently mis-joins: isec reports false discordance, merge emits duplicate rows and splits the allele frequency, dedup misses the duplicate. Nothing errors -- the counts are simply wrong.

Decision: **normalize (left-align + parsimony against the SAME reference, split multiallelics) every input before any merge/concat-dedup/isec/annotate.** This matters most for indels in homopolymers/STRs, where callers legitimately disagree on POS. It is safe to skip only when a single caller produced all inputs and no cross-file matching or database lookup follows. The canonical incantation is `bcftools norm -m-any -f ref.fa`; see bio-variant-normalization for the full pipeline, MNP atomization, and the vt-vs-bcftools discordance -- do NOT re-derive that here, cross-reference it.

## merge vs concat vs isec (choose by what differs)

| Operation | Inputs differ in | Produces | Requires index | Fails / misused when |
|-----------|------------------|----------|----------------|----------------------|
| `bcftools merge` | **samples** (same sites) | one multi-sample VCF (columns unioned) | yes | given the SAME sample split by region -> use concat; given single-sample VCFs and treated as joint genotyping -> fabricates 0/0 (see trap below) |
| `bcftools concat` | **regions** (same samples, e.g. per-chromosome) | one VCF spanning all regions (rows appended) | only with `-a` | inputs from DIFFERENT samples -> use merge; inputs overlap without `-a`; using `--naive` when headers/sample order differ |
| `bcftools isec` | neither -- same cohort, compare membership | per-input private/shared partition dirs | yes | inputs not normalized to identical representation -> false discordance |

Common confusion: `concat -a` (`--allow-overlaps`) resolves duplicate records from the SAME sample across overlapping region files; it does NOT union genotypes across different samples -- that is merge. If bcftools reports a "different samples" error, the operation is inverted.

## The merge trap: `bcftools merge` is not joint genotyping

When merging single-sample VCFs, a site called in sample A but simply **absent** from sample B's file is ambiguous: was B confidently homozygous reference there, or was B never covered/called? A project VCF cannot answer this. `bcftools merge` fills B's genotype with `./.` (missing) by default, and `-0/--missing-to-ref` overrides it to `0/0` -- but **both are guesses**, because merge has no evidence for the unseen site. `-0` therefore **fabricates** hom-ref genotypes and inflates the reference-allele count (see the vcf-basics `./.`-is-not-`0/0` distinction). Use `-0` only when every input truly covered every site (e.g. gVCF-derived, or a shared target with confirmed coverage), never as a convenience to remove `./.`.

Decision: to build a multi-sample callset with correct hom-ref-vs-no-data resolution, **joint-genotype gVCFs** (`GenomicsDBImport`/`CombineGVCFs` -> `GenotypeGVCFs`), do not `bcftools merge` single-sample project VCFs -- see upstream-only joint-calling. Reserve `bcftools merge` for combining already-jointly-genotyped cohorts, or samples that share a target with known coverage.

Merge also requires two harmonizations, or it silently drops or mis-collapses records:
- **Consistent representation.** All inputs must be normalized and split the same way first. If cohort A is split biallelic and cohort B keeps multiallelics, merge mis-collapses the shared site. Normalize all inputs identically (governing principle above).
- **Matching `##contig` headers and sample names.** Merge unions sample columns; duplicate sample names abort unless `--force-samples` renames them, and mismatched contig naming (`chr1` vs `1`) prevents sites from aligning. Rename samples with `bcftools reheader -s`; rename record CHROM values with `bcftools annotate --rename-chrs map.txt` (reheader only rewrites the `##contig` header lines, so a chr1-vs-1 merge after reheader silently doubles every site).

## bcftools merge (combine different samples)

```bash
# Union samples across per-sample (already joint-genotyped or shared-target) VCFs
bcftools merge -l files.txt -Oz -o cohort.vcf.gz    # -l: one VCF path per line
bcftools index cohort.vcf.gz
```

- `-m, --merge` controls multiallelic collapse at shared sites (default `both`): `-m none` keeps a SNP and an indel at one POS as separate records; `-m snps`/`-m indels` restrict which types collapse. Leave the default unless a downstream tool needs types kept apart.
- `--force-samples` disambiguates colliding sample names; `-r chr:beg-end` restricts to a region (inputs must be indexed).

## bcftools concat (stitch regions for the same samples)

```bash
# Genome-wide file from per-chromosome calls (same samples, disjoint regions)
bcftools concat chr{1..22}.vcf.gz chrX.vcf.gz -Oz -o genome.vcf.gz
```

- `-a, --allow-overlaps` is needed when region files overlap (e.g. windowed calling); pair with `-d/--rm-dups <snps|indels|both|all|exact>` to output a duplicate once. `-a` requires indexed inputs.
- `-n, --naive` concatenates BCF/VCF blocks WITHOUT recompression -- very fast for a large per-chromosome set, but it does only a header-compatibility check and requires identical headers and identical sample order across all files; it cannot reorder or reconcile anything. `--naive-force` skips even the header check and will silently produce a corrupt file if headers differ -- avoid it unless the files were provably produced identically.
- concat does NOT sort across file boundaries; overlapping unsorted inputs need `-a`, and the final file may still need `bcftools sort`.

## bcftools sort (order by CHROM then POS)

```bash
bcftools sort -T /scratch/tmp -m 4G input.vcf.gz -Oz -o sorted.vcf.gz   # -T tempdir, -m spill threshold for large files
```

An unsorted VCF breaks everything downstream: `tabix`/`bcftools index` require coordinate-sorted input to build the index, and merge/isec/`view -r` all depend on that index for random access. Sort after any operation that can leave records out of order (naive concat of misordered files, some reheader edits). `-T`/`-m` bound memory for genome-scale files.

## bcftools isec (set operations on call sets)

```bash
# Normalize BOTH first (governing principle), then partition
bcftools norm -m-any -f ref.fa gatk.vcf.gz     -Oz -o gatk.norm.vcf.gz
bcftools norm -m-any -f ref.fa freebayes.vcf.gz -Oz -o fb.norm.vcf.gz
bcftools isec -p comparison -Oz gatk.norm.vcf.gz fb.norm.vcf.gz
```

`-p dir` writes the four-way partition (`-Oz` to compress):

| File | Contents |
|------|----------|
| `0000.vcf[.gz]` | private to file 1 |
| `0001.vcf[.gz]` | private to file 2 |
| `0002.vcf[.gz]` | shared, file-1 records (file-1 INFO/FORMAT) |
| `0003.vcf[.gz]` | shared, file-2 records (file-2 INFO/FORMAT) |

`0002` and `0003` are the SAME sites with each file's own annotations -- pick by which annotations are needed downstream. Select membership instead of the full partition with `-n` and route records with `-w` (1-based file indices):

| Flag | Meaning |
|------|---------|
| `-n=2 -w1` | present in exactly 2 files, output file-1 records |
| `-n+2 -w1` | present in >=2 files |
| `-n~10 -w1` | present in file1 but NOT file2 (boolean mask) |
| `-C` | complement: positions only in file1, missing in the rest |

`-c, --collapse` sets what counts as "the same record"; the default `none` demands an exact REF+ALT match (why normalization is mandatory first), whereas `-c all` matches on position alone and ignores ALT -- rarely what a caller comparison wants.

## Subsetting samples and regions (`bcftools view`)

```bash
bcftools view -s sample1,sample2 input.vcf.gz -Oz -o subset.vcf.gz   # -s ^s3 to EXCLUDE; -S file for a list
bcftools view -r chr1:1e6-2e6      input.vcf.gz -Oz -o region.vcf.gz  # -R file.bed for many regions
```

Two nuances that bite:
- **`-r`/`-R` (regions) vs `-t`/`-T` (targets).** `-r`/`-R` use the index to JUMP to regions (fast, require an index) and consider both POS and an indel's end; `-t`/`-T` STREAM the whole file filtering on POS (no index needed, slower). Older bcftools releases could emit a record more than once for overlapping `-R` regions; 1.21 and 1.24 do not. If the version is unknown, merge overlapping intervals first.
- **Stale INFO counts after subsetting.** Dropping samples makes INFO `AC/AN/AF` wrong. `bcftools view -s` updates `AC/AN` by default (unless `-I/--no-update`), but recompute the full tag set explicitly: `bcftools +fill-tags subset.vcf.gz -Oz -o out.vcf.gz -- -t AC,AN,AF`.

## Header harmonization (`bcftools reheader`)

```bash
printf 'old_name\tnew_name\n' > rename.txt
bcftools reheader -s rename.txt input.vcf.gz -o renamed.vcf.gz   # -s renames samples only, no record rewrite
```

`reheader` rewrites only the header (fast, no record pass): `-s` maps sample names, `-h` swaps in a whole new header, `-f ref.fa.fai` fixes `##contig` lines (names and lengths) to match a reference but does NOT rename the CHROM column of records. For `chr1`-vs-`1` naming use `bcftools annotate --rename-chrs map.txt in.vcf.gz -Oz -o out.vcf.gz` (map: `old<TAB>new` per line). Harmonize sample names and record contig names BEFORE merge so columns and sites align.

## Structural variants merge differently -- do NOT use `bcftools merge`

For SVs (`<DEL>`/`<DUP>`/`<INV>`/BND), "the same event" is fuzzy: breakpoints disagree by CIPOS/CIEND margins, so tuple-exact bcftools operations treat one deletion called by two tools as two variants. SV merging needs coordinate-and-size (ideally sequence) aware tools -- Truvari, SURVIVOR, or Jasmine -- whose distance/size parameters ARE the result. Use bcftools here only for small variants; route SV consensus to upstream-only structural-variant-calling.


## References

- Danecek P, Bonfield JK, Liddle J, et al. Twelve years of SAMtools and BCFtools. *GigaScience.* 2021;10(2):giab008. doi:10.1093/gigascience/giab008 (bcftools merge/concat/isec/norm/view/reheader reference implementation)
- Tan A, Abecasis GR, Kang HM. Unified representation of genetic variants. *Bioinformatics.* 2015;31(13):2202-2204. doi:10.1093/bioinformatics/btv112 (why normalization before tuple-keyed merge/isec is mandatory)
