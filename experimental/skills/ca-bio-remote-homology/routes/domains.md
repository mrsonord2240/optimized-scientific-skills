# Profiles and domains

```bash
hmmscan --cut_ga --domtblout query.domtbl Pfam-A.hmm query.fa
```

Provide a locally pressed Pfam database and query. `--cut_ga` uses model gathering thresholds. The source `examples/pfam_annotation.sh query.fa pfam_dir` downloads and presses a database if missing and writes fixed output names; inspect storage and use a fresh directory before running. The source's organism-specific toy example and its bundled public fixtures are omitted from this generic draft; they remain in the source skill. For profile-profile comparison build a profile with `hhblits -i query.fa -d uniref30 -oa3m query.a3m -n 3 -cpu 8`, then `hhsearch -i query.a3m -d pdb70 -o query.hhr -cpu 8`; provide compatible databases and check versions. A known HMM can be searched against sequences using `hmmsearch profile.hmm database.fa`.

Done when model/database releases, thresholds, alignments and result statistics are recorded.
