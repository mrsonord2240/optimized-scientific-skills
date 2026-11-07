# Choose the assay-specific method

```bash
python scripts/cfdna_cli.py --help
```

This command displays the available consensus and QC operations. There is no turnkey minimal sWGS or native-end fragmentomics pipeline here: select and validate the trim, align and deduplication policy for the assay before proceeding. Record library chemistry, UMI design, input mass, sensitivity versus specificity, downstream goal, and any size selection. Simplex removes many amplification errors; duplex can suppress strand-specific damage but costs recovery and needs validation. Without UMIs avoid position-only deduplication for quantitative readouts. See [method selection](../references/method-selection-and-fragment-qc.md) for low-input and contamination caveats.

**Done when:** the chosen method and its assay-specific tradeoffs are documented, or the unsupported workflow is handed to an appropriate pipeline.
