# Background (not needed to run a route)

Every LFC, QC gate and hit call is computed against a reference committed once, at library-order time. A wrong but silent commitment invalidates the endpoint with no error thrown.

## Why the four commitments matter

| Commitment | What it feeds downstream |
|------------|--------------------------|
| Guide library (guide-gene map, non-targeting, CEGv2 essentials, NEGv1 non-essentials) | The counting denominator, the QC calibrator, the hit-calling priors. Non-targeting guides (about 1% of the library) calibrate the null and FDR; CEGv2 and NEGv1 calibrate PR-AUC and the BAGEL2 and Chronos priors. Swapping or dropping a class silently breaks one of them. |
| Baseline (plasmid pool, Day-0, vehicle) | Every LFC. Plasmid is the cloning-bottleneck baseline, Day-0 the biology baseline. Drug against Day-0 conflates drug effect with normal proliferation. |
| Screen type (dropout, enrichment, FACS, drug-modifier) | Which hit-calling method is valid at all. |
| Copy-number profile (cancer lines) | Whether amplicon artifacts leave before hit calling; residual correlation of LFC with CN is the tell. |

## Copy-number artifact

Several simultaneous Cas9 cuts at an amplified locus trigger a gene-independent DNA-damage arrest, so amplified regions look essential regardless of gene function. It appears in both TP53-mutant and wild-type lines, though its magnitude correlates with TP53 status; separately, p53-dependent toxicity of Cas9 cutting in general is reported. CRISPRcleanR corrects it unsupervised without a CN profile; Chronos takes CN as input; CRISPRi avoids the cut.

## Notes on thresholds

- CEGv2 PR-AUC above 0.7 is a community convention, not a threshold from the CEGv2 paper.
- Replicate Pearson 0.8 is the MAGeCK-VISPR floor. On real HAP1 TKOv3 data (T0 plus three T18 replicates), averaging every off-diagonal pair gave 0.677 and the true T18 replicate pairs 0.789, a gap straddling the floor.
- BAGEL2 BF above 6 corresponds to FDR below 3% in the Hart 2017 calibration; BF above 3 to FDR below 5%.
- BAGEL2 `bf` unseeded on HAP1 TKOv3: mean absolute BF difference 1.31, maximum 89.3 over 18,053 genes, 39 genes flipping the BF above 6 call. Rerun with `--seed 42` here: maximum difference 0.0.
- Single-hit calls at FDR 0.05 carry about 5% false discovery; tier-1 consensus shrinks it.
- Single-cell: MOI 0.3 gives about 26% infected cells and about 4% multi-infected.
- Batch is a covariate in MAGeCK MLE, not a ComBat pre-correction, because ComBat on counts distorts the negative-binomial mean-variance relationship the callers assume.
- Run-time note: `mageck test` on the 70,754-guide HAP1 table took about 1m45s.
