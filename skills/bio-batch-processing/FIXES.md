# Fix ledger

Audit source: `initial-opt10-20260928`, candidate
`244c82a53a96dce4586308b71ed678bae2bd3bed685a1db7aa44c6459b0c8990`.

| Finding | State | Implementing surface | Focused regression |
|---|---|---|---|
| BATCH-001 | fixed | `scripts/batch_process.py:split_by_prefix` validates prefixes, contains resolved targets, rejects cross-platform collisions, and creates outputs exclusively | traversal and existing-target tests |
| BATCH-002 | fixed | bounded LRU handle cache, unconditional cleanup, complete/partial manifest | 2,048-prefix bounded-handle test |
| BATCH-003 | fixed | explicit `SUMMARY_FIELDS` and header-only empty result | empty-summary CSV test |
| BATCH-004 | fixed | positive-integer validation before parser/output creation | zero, negative, float, Boolean, and one-record tests |
| BATCH-005 | fixed | importable worker plus guarded `main()` recipe | real two-worker spawn subprocess |
| BATCH-006 | fixed | `stable_paths` and stable sorting at every multi-file script boundary | opposite-creation-order byte comparison |
| BATCH-007 | fixed | routed `scripts/pyfastx_index.py` gzip-index example | build/reuse/random-access test |

Verification commands and exact outputs are recorded in
`F:\OpenScience\audits\bio-batch-processing\fix-opt10-20260928`.
