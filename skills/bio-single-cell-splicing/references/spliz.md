# SpliZ: annotation-free cell-state discovery

SpliZ computes a per-gene splicing Z-score and tests cell-state association. It is a Nextflow pipeline, not a Python package or standalone CLI.

```bash
nextflow run salzmanlab/spliz -r main -latest -c spliz.config
```

Configure `dataname`, `input_file`, `libraryType`, `grouping_level_1`, and `grouping_level_2`; choose SICILIAN output (`SICILIAN=true`) or BAM samplesheet, metadata, and GTF (`SICILIAN=false`). This command was not executed in the reference environment because Nextflow was unavailable. Verify current configuration keys in the checked-out pipeline before running.
