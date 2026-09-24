# Contig Naming, GRCh38 Flavours and Renaming

Moved from `SKILL.md` (2026-09-21); the rename-by-map recipe is now `scripts/rename_contigs.sh`. Read when a BAM and its reference disagree on contig names (`chr22` vs `22`, RefSeq accessions), when asked which GRCh38 a BAM was aligned to, or before renaming contigs without re-aligning.

### GRCh38 Is Not One Reference

Contig sets checked 2026-09-20 from the `.fai` / `chrom.sizes` files and NCBI's `README_analysis_sets.txt` (`ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/`).

| Reference flavor | ALT | Decoy | EBV | HLA | Use case |
|------------------|-----|-------|-----|-----|----------|
| UCSC `hg38.fa` (455 contigs) | yes (261 `_alt`) | no | no | no | UCSC browser tracks |
| NCBI `no_alt_analysis_set` | no | no | yes (`chrEBV`) | no | Aligners that are not ALT-aware |
| NCBI `no_alt_plus_hs38d1_analysis_set` | no | yes | yes | no | Same, with decoys |
| NCBI `full_analysis_set` | yes | no | yes | no | ALT-aware BWA-MEM |
| NCBI `full_plus_hs38d1_analysis_set` | yes | yes | yes | no | ALT-aware BWA-MEM with decoys |
| 1000G `GRCh38_full_analysis_set_plus_decoy_hla` = bwakit hs38DH; Broad `Homo_sapiens_assembly38.fasta` has the same 3,366 contigs | yes (261) | yes (2,385) | yes | yes (525 `HLA-*`) | GATK Best Practices, 1000G BAMs |
| T2T-CHM13 v2.0 | n/a | n/a | n/a | n/a | Distinct coordinates -- NOT interchangeable |

Mixing no-alt and ALT-aware BAMs in one cohort produces inconsistent multi-mapping behavior at HLA, KIR, and segmental-duplication regions. Standardize before joint calling.

### Contig Naming: The Silent Killer

| Convention | Source | chr1 | mitochondrion |
|-----------|--------|------|---------------|
| UCSC hg38 | UCSC Genome Browser | chr1 | chrM (16,569 bp, same sequence as Ensembl `MT`) |
| UCSC hg19 | UCSC Genome Browser | chr1 | chrM (16,571 bp, NC_001807 -- NOT the sequence of GRCh37 `MT`, 16,569 bp) |
| Ensembl (GRCh37, GRCh38) | Ensembl | 1 | MT |
| NCBI RefSeq FASTA (`GCF_*_genomic.fna`) | NCBI | NC_000001.11 (GRCh38), NC_000001.10 (GRCh37) | NC_012920.1 |
| NCBI analysis sets, 1000G GRCh38 analysis set, Broad hg38 | NCBI, 1000G, Broad | chr1 | chrM |
| 1000G phase 3 GRCh37 (`hs37d5`) | 1000G | 1 | MT (plus `GL*`, `NC_007605` EBV, `hs37d5` decoy) |

The RefSeq FASTA has no `chr1`; the `_assembly_report.txt` beside it maps every name (tab-separated, CRLF line ends): column 1 Sequence-Name (`1`, `X`, `MT` = Ensembl names for the chromosomes only; scaffolds get GRC names there, Ensembl uses their GenBank accession, e.g. `KI270706.1`), 5 GenBank (`CM000663.2`), 7 RefSeq (`NC_000001.11`), 10 UCSC-style-name (`chr1`). Use it as the rename map for `_alt`, `_random` and `chrUn` contigs; stripping `chr` only covers chr1-22, X, Y, M.

A BAM with `@SQ SN:chr1` cannot be analyzed against a `1`-named reference (and vice versa). Detect:
```bash
samtools view -H sample.bam | grep '^@SQ' | head -3
samtools dict ref.fa | head -3
```

### Rename Contigs Without Re-aligning

Renaming is a header-only change: BAM records store a contig index, so `samtools reheader` gives records identical to the input except the name. Only do it when the sequences are the same (compare `LN`, and `M5` where the BAM header has it -- a UCSC hg19 `chrM` is not GRCh37 `MT`).
```bash
# map.tsv: old<TAB>new, one contig per line (e.g. chr22<TAB>22, chrM<TAB>MT), or UCSC -> RefSeq from the assembly report:
grep -v '^#' GCF_000001405.40_GRCh38.p14_assembly_report.txt | tr -d '\r' | awk -F'\t' '$10!="na" && $7!="na"{print $10 "\t" $7}' > map.tsv   # a few UCSC contigs (chrUn_KI270752v1) have no RefSeq accession ("na"): reheader stops with "Duplicate entry na" if they stay
scripts/rename_contigs.sh sample.bam map.tsv renamed.bam   # header-only rename + index; a failed reheader leaves no renamed.bam
# add ref.fa as a 4th argument to also check names and lengths against ref.fa.fai (prints OK)

# UCSC -> Ensembl for the primary chromosomes only (hg38/GRCh38, not hg19); chr1_KI..._alt, chrUn_... etc. keep their names
samtools view -H sample.bam | sed -E -e 's/^(@SQ\tSN:)chrM\t/\1MT\t/' -e 's/^(@SQ\tSN:)chr([0-9]+|X|Y)\t/\1\2\t/' > renamed.hdr
samtools reheader renamed.hdr sample.bam > renamed.bam && samtools index renamed.bam
```
`SA:Z:`, `XA:Z:` and `OA:Z:` tags hold contig names as text and keep the old names; drop them with `samtools view -b -x SA -x XA renamed.bam` if a downstream tool reads them. `samtools reheader` also accepts a CRAM (checked on 1.24: exit 0, header `M5` kept, records decode identically with `-T` the renamed reference; decoding without `-T` needs a `REF_PATH`/cache entry for that `M5`). For VCF use `bcftools annotate --rename-chrs map.tsv`.
