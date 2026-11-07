# Resource files and clusters

```bash
python scripts/query.py uniref --cluster CLUSTER_ID --out cluster.json
```

Replace the cluster ID. For a public proteome, `python scripts/query.py proteome --id PROTEOME_ID --out proteome.fasta.gz` downloads the source client's compressed stream and overwrites the named output; estimate size and select a new path first. This query includes reviewed and unreviewed entries. The source distinguishes UniProtKB, UniRef, UniParc and Proteomes; it does not provide a dedicated UniParc client. See [API notes](../references/api.md).

Done when the requested cluster or file is checked for identity, completeness and release.
