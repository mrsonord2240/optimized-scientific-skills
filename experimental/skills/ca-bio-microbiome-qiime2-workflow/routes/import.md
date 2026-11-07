# Import paired reads

```bash
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest.tsv --input-format PairedEndFastqManifestPhred33V2 --output-path demux.qza
```

The V2 manifest is tab-separated with `sample-id`, `forward-absolute-filepath`, and `reverse-absolute-filepath`; verify actual encoding and paths. Then run `qiime demux summarize --i-data demux.qza --o-visualization demux.qzv` to inspect quality. Casava directory, EMP multiplexed reads and BIOM use different import types and formats; select and verify the installed release's help before adapting. EMP requires a demultiplexing step after import. Do not treat a successful import alone as evidence that quality encoding is correct.

Done when IDs, semantic type, encoding, sample counts and quality summary are checked.
