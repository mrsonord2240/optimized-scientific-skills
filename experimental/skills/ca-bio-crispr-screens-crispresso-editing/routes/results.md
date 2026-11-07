# Parse and compare results
~~~bash
python scripts/parse_crispresso.py sample_results/CRISPResso_on_sample_id > metrics.json
~~~
The parser reads mapping statistics and editing quantification; it computes aligned/input percentage and optionally includes metadata. Input read count must be nonzero. For comparison of existing runs:
~~~bash
CRISPRessoCompare crispresso_results/CRISPResso_on_control crispresso_results/CRISPResso_on_edited --output_folder crispresso_results/comparison
~~~
Inspect allele tables and alignment, and treat recurring deletions as possible rather than proven MMEJ. The source example scripts/crispresso_analysis.sh contains fixed example sequences; replace them before use. See [interpretation](../references/interpretation.md).
Done when metrics, comparability and allele-level limitations are reviewed.
