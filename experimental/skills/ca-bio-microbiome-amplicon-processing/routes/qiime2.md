# QIIME2 equivalents

No script. The QIIME2 path installs as its own conda environment (`qiime2-amplicon-<release>`).

```bash
qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed.qza \
    --p-trunc-len-f 220 --p-trunc-len-r 200 \
    --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza
```

- Variants: `denoise-single`; `denoise-pyro` for 454 and Ion Torrent; `denoise-ccs` for PacBio CCS (`--p-front`, `--p-adapter`, `--p-min-len`, `--p-max-len`).
- Trim primers in QIIME2 with `qiime cutadapt trim-paired --p-front-f FWD --p-front-r REV --p-discard-untrimmed`.
- Same rules as R: truncLen budget (`routes/truncation.md`), per-run denoising, merge tables afterwards.

Deblur (static positive filter, one fixed length, 16S only):

```bash
qiime deblur denoise-16S --i-demultiplexed-seqs qc.qza \
    --p-trim-length 250 --p-sample-stats \
    --o-representative-sequences rep-seqs.qza --o-table table.qza --o-stats stats.qza
```

`--p-trim-length` is one fixed value; Deblur cannot handle variable length. Do not merge a Deblur table with a DADA2 table.

Done when the table and representative sequences artifacts exist.
