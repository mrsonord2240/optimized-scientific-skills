# Calculate a TSS-relative window
~~~bash
python scripts/tss_windows.py --mode crispri --tss 1000000 --strand +
~~~
Supply a TSS coordinate, strand and mode. The source Dolcetto window is -50 to +300 relative to TSS; Calabrese CRISPRa is -150 to -75. Confirm coordinate convention and actual TSS from appropriate data; this computes bounds only and does not select or score guides. See [TSS notes](../references/crispri-crispra-tss.md).
Done when bounds, assembly, strand, TSS provenance and potential guide shortfall are recorded.
