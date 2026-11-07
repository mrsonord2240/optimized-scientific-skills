# Fungal ITS

ITS length is biological (ITS1 about 200-600 bp), so a fixed truncLen slices long variants and breaks merging of short ones. Trim with cutadapt, then ITSxpress, then denoise single-end with `--trunc-len 0`.

```bash
# primers first: routes/remove-primers.md
itsxpress --fastq r1.fastq.gz --fastq2 r2.fastq.gz \
    --region ITS2 --taxa Fungi \
    --outfile trimmed.fastq.gz --threads 4

Rscript scripts/dada2_asv.R --single --trunc-len 0 --min-len 50 --out its_out its_trimmed_dir
```

`--region` is ITS1, ITS2 or ALL; `--taxa` picks the HMM model. Paired input is merged, giving one trimmed FASTQ per sample. Put each run's ITSxpress outputs (`<sample>.fastq.gz`) in one directory per run. With `--single` the script filters on `maxEE` 2 and `minLen`, learns the error model per run, and does the same cross-run merge and single chimera removal.

Done when `its_out/seqtab_nochim.rds` exists.
