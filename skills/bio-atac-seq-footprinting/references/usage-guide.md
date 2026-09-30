# Footprinting prerequisites and request examples

The core Skill describes the workflow; [`method-reference.md`](method-reference.md) holds method details and caveats.

## Environments

The tools pin incompatible Pythons (rgt and pydnase need Python 3.7-era stacks, scPrinter needs 3.11 and pins), so use one environment per tool. Installing them together resolves TOBIAS 0.13.3 on Python 3.7, below the tested version, and `bioconda` without `conda-forge` cannot solve `tobias` (missing `adjusttext`). These are the commands the tested environments were built from (micromamba shown; conda or mamba work alike):

```bash
# TOBIAS 0.17.5 (also samtools, bedtools, deepTools alignmentSieve)
micromamba create -n footprint -c conda-forge -c bioconda python=3.10 pip pybigwig pysam samtools=1.19 bedtools numpy deeptools
micromamba run -n footprint pip install tobias

# HINT-ATAC (rgt 1.0.2) and Wellington (pyDNase 0.3.0), each alone
micromamba create -n footprint-rgt -c conda-forge -c bioconda rgt
micromamba create -n footprint-pydnase -c conda-forge -c bioconda pydnase

# scPrinter 1.2.0 (GitHub tag, not on PyPI): torch build for your CUDA, then pins
micromamba create -n footprint-scprinter -c conda-forge -c bioconda python=3.11 pip git bedtools
micromamba run -n footprint-scprinter pip install torch --index-url https://download.pytorch.org/whl/cu128   # tested: torch 2.11.0+cu128
micromamba run -n footprint-scprinter pip install "git+https://github.com/buenrostrolab/scPrinter@v1.2.0" "tangermeme==0.4.4" "snapatac2==2.8.0" ema_pytorch
```

The scPrinter pins exist because v1.2.0 imports `tangermeme.tools.tomtom` (absent in tangermeme 1.5.0) and `import_fragments` needs `snapatac2.pp.import_data` (absent in 2.9.0); snapatac2 2.8.0 downgrades numpy, pandas, and anndata and pip warns about zarr and shap, but the tested path ran. The pinned command was re-run in a fresh environment (about 30 min); it builds macs3, MOODS-python, and sorted_nearest from source, so it needs a C/C++ compiler and git. Set `SCPRINTER_DATA` to a writable directory; the first use downloads pretrained models and needs internet.

### RGT data for HINT-ATAC

The bioconda `rgt` package ships no data directory, and `rgt-hint` fails without `data.config`. Reinstall RGT from PyPI source with `RGTDATA` set once (the bioconda copy makes a plain `pip install` a no-op, so `--force-reinstall --no-cache-dir` is required; `setup.py` then writes `data.config` and copies the HMMs and bias tables), then point the hg38 entry at your own genome and annotation (or run `python setupGenomicData.py --hg38` to download them):

```bash
export RGTDATA=$HOME/rgtdata
pip install --no-deps --force-reinstall --no-cache-dir rgt==1.0.2   # in the footprint-rgt environment; writes $RGTDATA/data.config
cd "$RGTDATA" && python setupGenomicData.py --hg38 --hg38-genome-path hg38.fa --hg38-gtf-path gencode.gtf
```

### Fragment file for scPrinter

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

## Motif database

```bash
wget https://jaspar.genereg.net/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
mv JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt JASPAR2024_CORE_vertebrates.pfm   # default motif name in scripts/run_tobias.sh
```

## Request examples

- "Run the TOBIAS three-step pipeline (ATACorrect, ScoreBigwig, BINDetect) on two conditions, with `--cond-names treated control`."
- "Plot the aggregate footprint at JASPAR CTCF MA0139.2 sites with TOBIAS PlotAggregate and confirm a clean dip at bound sites."
- "Run scPrinter for multi-scale footprints on this bulk library and compare CTCF bound and unbound sites."
- "Score CTCF footprints per leiden cluster of this scATAC fragment file with scPrinter and report whether bound sites separate from unbound ones in each cluster."
- "Run TOBIAS and HINT-ATAC; report the bound-versus-unbound site overlap of the HINT footprints as the concordance measure."
- "Filter to fragments under 100 bp, run TOBIAS on the NFR BAM, and compare footprint sharpness with the full-fragment BAM."
- "Cross-validate predicted FOXA1 footprints against published ChIP-seq peaks."

Provide the BAMs per condition, consensus peaks, genome FASTA and build, blacklist, motif database, organism, and whether the data are bulk or single-cell. Ask for clarification when these choices change the tool.
