# Somatic variant filtering (Mutect2)

Step 1 needs GATK and the contamination and segmentation tables, so it has no wrapper:

```bash
gatk FilterMutectCalls -R reference.fa -V mutect2_raw.vcf \
    --contamination-table contamination.table \
    --tumor-segmentation segments.table \
    -O mutect2_filtered.vcf
```

Step 2, extra thresholds on the tumor column:

```bash
bash scripts/somatic_filter.sh mutect2_filtered.vcf somatic_final.vcf tlod=6.3 af=0.05 dp=20
```

- The script finds the tumor column from the `##tumor_sample` header (or `tumor=<name>`). Mutect2 often lists the normal first, so a fixed `[0]` would test the normal. It exits 2 when the tumor sample is not found, instead of letting bcftools segfault.
- The expression is a site filter on the tumor column; the normal column is not tested.
- By default it does not require `FILTER=PASS`, so records FilterMutectCalls rejected are kept if they meet the thresholds. Pass `pass_only=1` to drop them. Say which you used.
- Never use VQSR or the germline hard filters on somatic calls.

Done when the script printed the tumor sample name and record counts in and out.
