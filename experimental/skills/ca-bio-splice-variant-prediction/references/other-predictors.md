## SpliceTransformer, CI-SpliceAI and CADD-Splice

**SpliceTransformer** (You 2024, GitHub `ShenLab-Genomics/SpliceTransformer`; SNVs only, indels skipped silently; hg38 or hg19). Weights are not in the repo: a public Google Drive file (`gdown 1d8n4vHDSbXqpPc_JFEswLomSUDBgHvno`, 120 MB) saved as `model/weights/SpTransformer_pytorch.ckpt`.

```bash
pip install torch sinkhorn-transformer "axial-positional-embedding==0.2.1" "PyVCF3==1.0.0" pyensembl gffutils pyfaidx pandas tqdm gdown
# axial-positional-embedding 0.2.1 is required: newer versions rename pos_emb weights and load_state_dict fails
# put a chr-prefixed hg38.fa and GENCODE hg38.annotation.gtf.gz in data/data_package/, then index the GTF once
# (the annotation name is hard-coded in sptransformer.py, whatever GENCODE release the GTF is):
pyensembl install --reference-name hg38 --annotation-name gencode.v38 --gtf data/data_package/hg38.annotation.gtf.gz
python sptransformer.py -I input.vcf -O out.csv --reference hg38
```

Output: `score` (largest change in splice-site probability, 0-1) plus Y/N flags for 15 tissues (tissue-specific usage change above the paper's 95% threshold). Panel result: TP53 c.673-2A>G 1.00, DMD c.9563+1G>A 0.99, GLA c.370-1G>A 0.98, GLA c.639+919G>A 0.38, ClinVar-benign GLA 0.01-0.07. Runtime is ~25 s per variant on CPU (3 min 52 s for 9 variants) and 12 s for the same 9 on a GPU: the RTX 5070 Ti (Blackwell, sm_120) needs a CUDA 12.8 torch build (`pip install torch --index-url https://download.pytorch.org/whl/cu128`, checked with torch 2.11.0+cu128); scores were identical on CPU and GPU.

**CI-SpliceAI** (Strauch 2022; SpliceAI retrained on all GENCODE isoforms). `pip install "tensorflow-cpu==2.15.*" "keras<3" cispliceai` (keras 3 is not supported; 1.2.2 checked), then:

```bash
cis-vcf -a grch38 -d 500 --all -i input.vcf -o ci_output.vcf GRCh38.primary_assembly.genome.fa
```

Output INFO tag is `CISpliceAI=ALLELE|ENSEMBL_GENE|DS_AG|DS_AL|DS_DG|DS_DL|DP_...` (SpliceAI layout; parse it by replacing `SpliceAI=`); `-a grch37` for GRCh37. On the GRCh37 panel it flagged every canonical variant (delta_max 0.64-1.00; SpliceAI 0.95-1.00) but scored the GLA c.639+919G>A pseudoexon 0.05 at `-d 500` where SpliceAI gave 0.30, so it is a cross-check, not a rescue.

**CADD-Splice** is the standard CADD score from v1.6 on (splice predictors are inputs to the model). The public API needs no account:

```bash
curl -s "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.7/17:7674292_T_C"
# [{"Alt":"C","Chrom":"17","PHRED":"34","Pos":"7674292","RawScore":"6.512109","Ref":"T"}]
```

PHRED on the panel (GRCh37-v1.7): canonical splice variants 33-35 (DMD c.31+1G>A 33, DMD c.9563+1G>A 34, GLA c.370-1G>A 35), GLA c.639+919G>A 14.9, benign GLA rs2071228 10.7. It does not say which component drove the score; report it beside, not instead of, SpliceAI/MMSplice.
