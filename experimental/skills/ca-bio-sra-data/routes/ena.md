# Public ENA retrieval

```bash
bash examples/download_batch.sh runs.txt ./new-output-directory
```

Supply one known-public uppercase SRR, ERR or DRR run ID per line. Before retrieval estimate size and storage; source requests explicit confirmation before transfers above approximately 10 GB. The script queries ENA fields by header, downloads files and compares MD5s. It writes `failed.txt` and can exit zero despite failures; inspect that file. Use a new empty output directory: this batch script can overwrite existing names. Missing URLs or authorization errors require investigation, not automatic retries around access restrictions. Compare counts to ENA metadata, not another archive. See [access notes](../references/access.md).

Done when all intended public runs have verified files and no unresolved entries in `failed.txt`.
