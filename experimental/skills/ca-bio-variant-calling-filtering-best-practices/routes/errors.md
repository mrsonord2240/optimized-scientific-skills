# Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `no such INFO tag` | Tag absent from the VCF | Check the header: `bcftools view -h in.vcf` |
| `syntax error` in an expression | Invalid operator | Use `\|\|`, not `or`; write missing as `= "."` |
| Every hom-alt site removed | A RankSum term has no missing-value guard | Add `\|\| INFO/X = "."` to each RankSum term |
| VQSR "converged" but the scores are nonsense | Too few variants | Hard filters, VETS or NVScoreVariants |
| Empty output | Filter too strict | Relax thresholds; inspect annotation histograms |
| `chromosome block chr1 is not contiguous` | `concat` of position-interleaved SNP and indel files | Index both and use `bcftools concat -a` |
| `tag AF can have multiple subfields, run as AF[sample:subfield]` | A Number=A/R FORMAT tag indexed as `AF[0]` | Use `FMT/AF[T:0]` (T is the sample column) or `FMT/AD[:1]` for every sample |
| bcftools exit 139 on a somatic filter | Tumor column index -1 (no `##tumor_sample` header) | Use `scripts/somatic_filter.sh ... tumor=<name>` |

Done when the failing command has run to completion.
