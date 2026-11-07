# Callset statistics

```bash
bash scripts/vcf_qc.sh input.vcf.gz qc_out
PLOT=1 bash scripts/vcf_qc.sh input.vcf.gz qc_out      # also plot-vcfstats
```

Writes `stats.txt` (`bcftools stats`), `per_sample.txt` (`-s -`), `summary.txt` (SN counts, TSTV, per-sample PSC) and `counts.txt` (records, SNPs, indels, mean QUAL, PASS counts). Read the numbers with `routes/interpret.md`.

- Compare two callsets: `bcftools stats file1.vcf.gz file2.vcf.gz > cmp.txt`.
- `stats.txt` sections: `SN` summary, `TSTV`, `SiS` singletons, `AF` spectrum, `QUAL`, `IDD` indel lengths, `ST` substitution types, `DP`, `PSC` per-sample counts (hom-ref, het, hom-alt, transitions, transversions, missing), `PSI` per-sample indels.
- `FILTER=PASS` strictly excludes records with `FILTER='.'`; on an unfiltered callset use `bcftools view -f .,PASS`. `counts.txt` prints both.
- `plot-vcfstats` needs matplotlib for PNGs and pdflatex or tectonic for `summary.pdf`; without LaTeX it exits 2 at the PDF step but the PNGs are written.

## Python (cyvcf2)

```bash
python scripts/vcf_stats.py input.vcf.gz [--per-sample] [--afs]
```

Prints counts, FILTER-null versus filtered, allele-level Ti/Tv (a multiallelic SNP counts once per ALT) and mean QUAL; `--per-sample` adds het/hom-alt/missing per sample, `--afs` an `INFO/AF` spectrum.

Done when `summary.txt` (or the Python output) exists and you have read Ti/Tv, the per-sample het/hom and the missing counts.
