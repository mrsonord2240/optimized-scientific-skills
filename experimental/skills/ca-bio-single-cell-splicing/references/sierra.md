# Sierra: alternative polyadenylation in 10X 3' data

Sierra detects alternative polyadenylation (APA), not cassette-exon splicing. Supply a comprehensive GTF with 3' UTR annotation and a junction file.

```bash
regtools junctions extract -a 8 -m 50 -M 500000 -s RF -o junctions.bed possorted_genome_bam.bam
```

```r
library(Sierra)
FindPeaks(output.file='peaks.txt', gtf.file='annotation.gtf',
          bamfile='possorted_genome_bam.bam', junctions.file='junctions.bed')
CountPeaks(peak.sites.file='peaks.txt', gtf.file='annotation.gtf',
           bamfile='possorted_genome_bam.bam', whitelist.file='barcodes.tsv', output.dir='peak_counts/')
counts <- ReadPeakCounts(data.dir='peak_counts/')
AnnotatePeaksFromGTF(peak.sites.file='peaks.txt', gtf.file='annotation.gtf',
                     output.file='peak_annotations.txt')
annot <- read.table('peak_annotations.txt', header=TRUE, sep='\t', row.names=1)
peaks.sce <- NewPeakSCE(peak.data=counts, annot.info=annot, cell.idents=cell_identities,
                        min.cells=0, min.peaks=0)
apa_results <- DUTest(peaks.sce, population.1='ctrl', population.2='trt')
```

`CountPeaks` and `AnnotatePeaksFromGTF` write files and return `NULL`. With SeuratObject >=5, use `NewPeakSCE`; `NewPeakSeurat` can fail because the `slot` argument is defunct. When donors exist, pass donor-defined `replicates.1` and `replicates.2` to `DUTest` rather than relying on random splits. Missing 3' UTR annotation means missed peaks, not evidence of absent APA.
