# Germline evidence ledger
~~~sh
python scripts/evidence_ledger.py evidence.json result.json
~~~
Use only synthetic or public non-patient examples. The input JSON has a criteria array, for example {"criteria":["PM2_Supporting"]}; entries are pre-reviewed evidence codes, not automatic annotations. The script checks source-implemented code conflicts, applies PVS1 subsumption and returns a provisional point-based category with a non-diagnostic notice. Preserve rationale, excluded evidence, conflict review, disease mechanism, transcript and the exact rule version separately; this helper cannot validate those inputs. Check the current gene/disease-specific VCEP specification and ClinGen SVI guidance before any qualified review. The [framework](../references/acmg-framework.md) documents the source thresholds and limitations; they are historical and not independently verified here.
**Done when:** output and unresolved context are retained for qualified review, with no autonomous clinical use.
