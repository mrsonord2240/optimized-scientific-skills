# Scprinter



Raw, unshifted fragments of properly paired reads, coordinate-sorted, bgzipped and indexed with the `footprint` environment's samtools and htslib. scPrinter's shift parameters give the shift already applied to the ends, so these take `--shift 0,0` (the script default; scPrinter's auto-detection also returned 0,0 on them) and Cell Ranger fragments, already shifted, take `--shift 4,-5`:

```bash
samtools view -f 2 sample.bam | awk 'BEGIN{OFS="\t"} $9>0 && $9<1000 {print $3,$4-1,$4-1+$9,"sample"}' \
    | sort -k1,1 -k2,2n -k3,3n > frags.tsv && bgzip frags.tsv && tabix -p bed frags.tsv.gz
```

### scPrinter per-cluster footprints (scATAC)

`scripts/scprinter_footprint.py --groups` scores each cell group as one pseudobulk. Make the barcode-to-cluster TSV in a separate process from scPrinter, or import scprinter before snapatac2: the reverse order breaks `import wandb`.

```bash
# cells.tsv: barcode<TAB>cluster, no header (for example snapatac2 leiden labels)
micromamba run -n footprint-scprinter python scripts/scprinter_footprint.py --fragments fragments.tsv.gz \
    --fasta hg38.fa --gtf genes.gtf --blacklist hg38-blacklist.v2.bed --regions ctcf_sites.bed --outdir sc_out \
    --groups cells.tsv --shift 4,-5 --min-fragments 1 --modes 10,20,30,50
```

Tested on 10x PBMC 5k fragments, chr1:1-30 Mb: 697 cells in 5 snapatac2 leiden clusters (49-198 cells, about 650 fragments per cell) at 400 CTCF motif sites labelled bound or unbound by TOBIAS. With an existing `--bias` file the run took 2.5 min. Scores were finite, and the mode-10 profiles at bound sites differed between clusters (pairwise r 0.04-0.59). The centre score at bound sites exceeded unbound in all 20 cluster-by-mode cells, but only 1 reached p < 0.05 (one-sided Mann-Whitney). At this depth the route runs but per-cluster footprints are not resolved; pooling all 2,973 barcodes with at least 200 fragments reached p 0.012 at mode 10 in the tooling run. Do not use `--shift auto` on sparse 10x fragments: it detected an implausible 22/12 there.

### seq2PRINT model training (GPU)

The bulk seq2PRINT route trains a sequence model per sample. Its train, valid and test splits are whole chromosomes, and `scp.pp.call_peaks` runs an executable named `macs2`; the environment ships `macs3`, so link it first:

```bash
mkdir -p shim && ln -sf "$(micromamba run -n footprint-scprinter which macs3)" shim/macs2 && export PATH=$PATH:$PWD/shim
export SCPRINTER_DATA=$HOME/scprinter_data
micromamba run -n footprint-scprinter python prep_seq2print.py
```

`prep_seq2print.py` (raw fragments from the recipe above; `Genome(name=...)` must be the pickle path, otherwise the final attribution step fails with `genome not supported`):

```python
import os, pickle
import pandas as pd
import scprinter as scp
work, fa, gtf, bl, frags = os.path.abspath("s2p"), "hg38.fa", "genes.gtf", "hg38-blacklist.v2.bed", "frags.tsv.gz"
splits = [{"train": ["chr1", "chr3", "chr5"], "valid": ["chr2"], "test": ["chr4"]}] * 5
os.makedirs(work, exist_ok=True)
bias = f"{work}/bias.h5"
if not os.path.exists(bias):
    scp.genome.predict_genome_tn5_bias(fa_file=fa, save_name=bias, tn5_model=scp.datasets.pretrained_Tn5_bias_model,
                                       context_radius=50, device="cuda:0", batch_size=5000)
genome = scp.genome.Genome(name=f"{work}/genome.pkl", fa_file=fa, gff_file=gtf, bias_file=bias, blacklist_file=bl, splits=splits)
pickle.dump(genome, open(f"{work}/genome.pkl", "wb"))
printer = scp.pp.import_fragments(path_to_frags=frags, barcodes=None, savename=f"{work}/printer.h5ad", genome=genome,
                                  sample_names=["sample"], min_num_fragments=1000, min_tsse=0, auto_detect_shift=False,
                                  plus_shift=0, minus_shift=0, sorted_by_barcode=False, low_memory=False)
scp.pp.call_peaks(printer=printer, frag_file=frags, cell_grouping=[None], group_names=["all"], preset="seq2PRINT", n_jobs=1)
pd.DataFrame(printer.uns["peak_calling"]["all_cleaned"][:]).to_csv(f"{work}/peaks.bed", sep="\t", header=False, index=False)
scp.tl.seq_model_config(printer, region_path=f"{work}/peaks.bed", cell_grouping=list(printer.obs_names), group_names="sample",
                        genome=genome, fold=0, overwrite_bigwig=True, model_name="sample", path_swap=(work, ""),
                        additional_config={"epochs": 3, "early_stopping": 3, "batch_size": 64,
                                           "dispersion_model": str(scp.datasets.pretrained_dispersion_model_v2)},
                        config_save_path=f"{work}/cfg_fold0.json")
printer.close()
```

Train fold 0 (wandb is optional; without `--enable_wandb` it is not used). Training ends by saving the model and running the count and footprint DeepSHAP attribution normalisation:

```bash
WANDB_MODE=disabled micromamba run -n footprint-scprinter seq2print_train --config "$PWD/s2p/cfg_fold0.json" \
    --data_dir "$PWD/s2p" --temp_dir "$PWD/s2p/temp" --model_dir "$PWD/s2p/model"
```

Tested only as a bounded route check. The staged real data covered chr1:1-30 Mb, so it was cut into three pseudo-chromosomes (16, 7 and 7 Mb) for train, valid and test; that layout is a test device, not an analysis design. 521,478 GM12878 fragments. The cleaned seq2PRINT peak set (12,579 peaks) validated at about 30 s per batch, too slow for the 45 min cap, so the config's `peaks` entry was pointed at 1,500/250/250 train/valid/test peaks. With 3 epochs, batch 64 and the default 11.8M-parameter model: training plus attribution finished with rc 0 in 31 min on a 16 GB RTX 5070 Ti. It saved the model, which has finite parameters and a forward output of 99 scales x 800 bp, and wrote count and footprint DeepSHAP bigwigs. Validation profile Pearson stayed at or below 0.07, so the model learned little: the run shows the commands work, not that the model is usable. GPU use peaked at 15.7 of 16 GB, with about 2-3 GB held by other processes, so plan for 11-16 GB of free GPU memory (the tooling run measured the same range).

LoRA single-cell seq2PRINT (`seq_lora_model_config`) and `seq_tfbs_seq2print` were not run. They need a bulk model trained per fold on genome-wide data, multi-million-fragment pseudobulks and cell embeddings, which this test data could not supply.
