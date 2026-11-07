# Sample swaps, relatedness, sex

```bash
OUT=id_qc bash scripts/identity_qc.sh peddy input.vcf.gz cohort.ped
OUT=id_qc bash scripts/identity_qc.sh somalier sites.vcf.gz ref.fa cohort.ped input.vcf.gz
OUT=id_qc bash scripts/identity_qc.sh king input.vcf.gz
OUT=id_qc bash scripts/identity_qc.sh gtcheck input.vcf.gz [truth.vcf.gz]
```

Build the all-pairs relatedness matrix and reconcile it with the manifest: a swap shows as unexpected relatedness, a failed "replicate" as a pair that is not related.

| Tool | Input | Detects |
|------|-------|---------|
| `gtcheck` | VCF | Same-file swaps and duplicates; a discordance score, no relatedness graph. Without a truth VCF it cross-checks all samples in the file |
| KING (`vcftools --relatedness2`) | VCF | Kinship without allele-frequency or ancestry assumptions |
| peddy | VCF + PED | Reported versus inferred sex, relationships, ancestry. Human genome-wide or exome only |
| somalier | BAM, CRAM, VCF | Relatedness, ancestry, sex at scale; cross-checks RNA-seq against WGS |

KING kinship bands: >0.354 duplicate or MZ twin; 0.177-0.354 first degree; 0.0884-0.177 second; 0.0442-0.0884 third. A pair not expected to be related at ~0.5 is a swap or duplicate.

- `bcftools gtcheck -G` no longer exists (removed in the 1.10 rewrite); the error probability flag is `-E`.
- Ancestry projection with somalier: `somalier ancestry --labels 1kg-labels.tsv 1kg/*.somalier ++ extracted/*.somalier`.
- peddy on a non-human or small targeted panel crashes or gives nonsense; use somalier or `gtcheck`.
- somalier `--sites` must match the genome build; build a custom one with `somalier find-sites <population.vcf.gz>`.

Done when the relatedness output exists and every unexpected pair or sex mismatch is listed.
