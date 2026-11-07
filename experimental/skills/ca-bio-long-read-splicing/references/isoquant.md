## IsoQuant for Discovery + Quantification

**Goal:** De novo or annotation-guided isoform discovery and quantification with high precision.

**Approach:** Run `isoquant` (not `isoquant.py`) with reference + reads (or a minimap2 BAM via `--bam`) + data type; output is GTF + counts.

```bash
isoquant \
    --reference reference.fa \
    --genedb gencode.v45.annotation.gtf \
    --fastq sample1.fastq.gz sample2.fastq.gz \
    --data_type pacbio_ccs \
    --output isoquant_output \
    --threads 16 \
    --model_construction_strategy default_pacbio
```

`--data_type` accepts `pacbio_ccs` (HiFi), `pacbio`, `nanopore`/`ont`, `assembly` or `transcripts`; `--stranded forward` for reads in transcript orientation (direct RNA), default `none`. `--genedb` is optional for de novo discovery. Outputs are in `<output>/<prefix>/`: `<prefix>.transcript_models.gtf` (also holds gene, CDS, UTR and codon rows copied from the annotation), `<prefix>.transcript_counts.tsv` (annotated transcripts only, pooled over all input files; default prefix `OUT`) and `<prefix>.discovered_transcript_counts.tsv` (annotated transcripts plus the novel models, ids like `transcript115.chrQ.nnic`), each with a `..._grouped_file_name_counts.tsv` twin (one column per input file; its id column is headed `gene_id` although the rows are transcripts). 4.0.0 writes no `transcript_model_counts.tsv`. Transcript counts default to `--transcript_quantification unique_only`, so ambiguous reads are not counted (`__ambiguous` row: 22 of 300 reads in the ONT microexon test). A novel splice site shifted by only ~24 nt from an annotated one is merged into the annotated isoform (all 20 planted reads of one sample counted in it, +13%); confirm such sites with short-read junctions (FLAIR `--junction_tab`). A planted 30-nt donor shift had no model in a `--genedb` run either, and `--illumina_bam` (125,055 simulated short reads, 40,278 spliced) changed nothing; a run without `--genedb` modelled it (43 reads), as it did an unannotated exon skip (37 reads, also modelled with `--genedb`). IsoQuant (Prjibelski 2023 *Nat Biotech*) is current SOTA for novel transcript reconstruction; pairs well with SQANTI3 for downstream classification.

Memory requirement: >=64 GB for atlas-scale runs.
