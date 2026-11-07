# SRA Toolkit and cloud retrieval

```bash
bash examples/download_single.sh RUN_ACCESSION ./new-output-directory 8 100G
```

Replace the placeholder with a known-public run. The source script validates accession syntax, uses prefetch with explicit size limit, validates SRA data, extracts and compresses files, and may fall back to ENA; review storage needs and the access boundary first. It skips technical reads. The alternative `bash examples/prefetch_large.sh RUN_ACCESSION ./new-output-directory 8 yes` retains technical reads and first tries the AWS public mirror; check cloud region and charges before use. Both clean their owned staging directories and refuse existing published names. Confirm size before transfers above approximately 10 GB; uncompressed extraction needs substantial scratch. Verify toolkit versions and flags, and inspect source-specific integrity and output. See [access notes](../references/access.md).

Done when the requested public run is validated and the expected files are present, or the failure is explicitly recorded.
