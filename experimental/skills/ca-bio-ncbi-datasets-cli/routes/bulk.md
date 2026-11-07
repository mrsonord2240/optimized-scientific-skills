# Bulk genomes

```bash
datasets download genome taxon "Bacteria" --reference --annotated --assembly-source RefSeq --include genome,gff3,protein --dehydrated --filename bact_refs.zip --no-progressbar
```

Inspect the ZIP and its `ncbi_dataset/fetch.txt`, unzip into a new directory, then run `datasets rehydrate --directory bact_refs --max-workers 8 --no-progressbar`. The source also provides `examples/bulk_dehydrated.sh Bacteria bact_refs.zip bact_refs 8` for aria2c transfer and size-based retry; review it before running: it deletes files whose sizes differ from the catalog. Verify fetch column layout and catalog schema for the installed version. Rehydrate may skip existing truncated files; size checks do not establish content integrity. See [bulk details](../references/dehydrated-bulk.md).

Done when requested assemblies are present and transfer integrity has been independently checked.
