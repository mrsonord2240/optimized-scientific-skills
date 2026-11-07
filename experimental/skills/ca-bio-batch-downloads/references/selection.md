# Selection



| Record count | Source | Strategy | Why |
|---|---|---|---|
| < 200 known IDs | Any db | EFetch with comma-joined `id=` | Single round-trip; trivial |
| 200-5,000 known IDs | Any db | EPost (chunked at 200) -> history -> chunked EFetch | URL length limit + chunked retrieval |
| 5,000-100,000 from a query | Any db | ESearch with `usehistory='y'` -> chunked EFetch | Push to server once; pull in batches |
| > 100,000 sequences | nucleotide/protein | Consider FTP mirror or Datasets CLI; if E-utils is still the chosen path, chunk via the history server (single stream) | NCBI throttles bulk; offline mirror is faster |
| > 1,000,000 sequences, or a literal "entire database" request | nucleotide/protein | **Refuse the EFetch-loop approach outright.** Point only to NCBI's bulk FTP/BLAST-db mirrors (`ftp.ncbi.nlm.nih.gov`) or Datasets CLI; do not attempt a chunked E-utilities loop at this scale | At this scale a single-stream E-utilities loop is a multi-week job this Skill's tools were never designed to serve -- not a "consider," a hard stop, the same way >4 concurrent workers is a hard stop |
| Whole genome assemblies | Assembly/Datasets | `datasets download genome accession ...` | Datasets v2 is the modern bulk endpoint |
| All RefSeq for a species | Datasets | `datasets download genome taxon ...` | Replaces assembly_summary.txt scraping |
| All gene records for a list | Datasets | `datasets download gene gene-id ...` | Cleaner output than EFetch gene XML |
| Raw sequencing reads | SRA | `prefetch` + `fasterq-dump` (or ENA mirror) | See `sra-data` skill |

The Datasets CLI is the right answer for any genome- or gene-centric bulk workflow as of 2023+. The E-utilities remain right for PubMed, ESummary metadata, custom queries, and anything not in the Datasets API. See `bio-ncbi-datasets-cli` skill.

Gene-related field tags (`[GENE]`, `[Gene Name]`) require the database's official gene symbol, not a descriptive term -- e.g. `insulin[Gene Name]` returns Count=0, `INS[GENE]` does not. Verify the symbol via ESearch or the `bio-entrez-search` skill's EInfo pattern before building a batch query around it.
