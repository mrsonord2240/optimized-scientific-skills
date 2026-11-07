# Batch, pooled and WGS
~~~bash
CRISPRessoBatch --batch_settings batch_settings.txt --batch_output_folder batch_run --skip_failed --n_processes 8
~~~
Prepare the tab-separated per-sample settings as in [multi-mode reference](../references/batch-pooled-wgs.md); inspect skipped failures. For pooled amplicons use:
~~~bash
CRISPRessoPooled --fastq_r1 pooled_R1.fastq.gz --fastq_r2 pooled_R2.fastq.gz --amplicons_file amplicons.txt --output_folder pooled_run --min_reads_to_use_region 100 --n_processes 8
~~~
For WGS use:
~~~bash
CRISPRessoWGS --bam_file aligned.bam --reference_file genome.fa --region_file off_targets.txt --output_folder wgs_run --n_processes 8
~~~
Inspect NA rows: the source reports default minimum 1000 reads for pooled and 10 for WGS, with skipped regions not necessarily failing the run. Verify table schemas, depth and installed version.
Done when each sample or region is accounted for and failures and NA results reviewed.
