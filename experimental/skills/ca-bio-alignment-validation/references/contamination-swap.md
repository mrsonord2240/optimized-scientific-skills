### Contamination and Sample Swap

No alignment QC is complete without these in production. `verifybamid2` and `somalier` need whole-genome or exome-scale data covering their marker panels (a region slice leaves too few markers: VerifyBamID2 2.0.3 reported "Insufficient Available markers" on a 100 kb slice).
```bash
# Cross-sample contamination; the prefix of the released SVD files includes ".dat"
verifybamid2 --SVDPrefix /resources/1000g.phase3.10k.b38.vcf.gz.dat \
    --Reference ref.fa --BamFile sample.bam --Output sample.contam
# FREEMIX > 0.03 is the commonly used contamination-concern threshold; values escalate from there.

# Relatedness, sex check, sample swap detection
somalier extract -d extracted/ -s /resources/sites.hg38.vcf.gz \
    -f ref.fa sample.bam
somalier relate --infer extracted/*.somalier

# Tumor/normal pairing verification: compare whole files so duplicate RG IDs/PUs cannot collapse
# two samples into one group. Alternatively give every BAM a distinct RG ID and PU.
# Give paired tumor/normal BAMs the same SM (patient id) or add EXPECT_ALL_GROUPS_TO_MATCH=true.
picard CrosscheckFingerprints I=tumor.bam I=normal.bam \
    HAPLOTYPE_MAP=Homo_sapiens_assembly38.haplotype_database.txt \
    CROSSCHECK_BY=FILE LOD_THRESHOLD=-5 OUTPUT=crosscheck.metrics
# With LOD_THRESHOLD=-5: LOD > 5 = same individual; < -5 = different; in between = ambiguous
```

Checked here: `CrosscheckFingerprints` and `somalier extract` + `relate --infer` ran end-to-end on a 40-site chr22 haplotype map / sites VCF built for the test slice (same reads: LOD +12.6, somalier relatedness 1.0; genotype-flipped reads: LOD -63, exit 1 for a same-SM pair, somalier relatedness -140, concordance 0). `verifybamid2` ran to the marker check only, so no FREEMIX value was produced here; the `HAPLOTYPE_MAP` above is the GATK-bundle file and the somalier/VerifyBamID2 panels are downloads, none shipped with this Skill.
