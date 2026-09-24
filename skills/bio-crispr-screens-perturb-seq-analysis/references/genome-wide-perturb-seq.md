## Genome-Wide Perturb-Seq (Replogle 2022)

**Replogle 2022 *Cell* 185:2559** demonstrated genome-wide Perturb-seq:
- >2.5M cells total; the genome-scale K562 screen targeted ~9,866 expressed genes (with a 2,057-gene essential subset)
- CRISPRi via dCas9-KRAB
- Native 10X 3' direct-capture for sgRNA
- Median >100 cells per perturbation as screened (Replogle 2022)
- Cluster-based analysis of perturbed cells reveals gene-program organization

**Scaling principles (calibrated to Replogle 2022's actual scope -- ~9,866 expressed genes -- not a generic "genome-scale" constant):**
- Cells per perturbation: 500-1,000 minimum for stable DE (Replogle's own screen ran leaner, at a median >100 cells/pert, trading DE power for genome-wide breadth -- see Failure Modes below)
- 10X channels: 10-30 channels at 5,000-10,000 cells each
- Cost: ~$50-100K

**Scaling to a different target gene count:** the figures above assume ~9,866 genes and scale roughly linearly: `figure x (target_genes / 9866)`. A literal ~19,000-protein-coding-gene design (this Skill's own usage-guide.md example) is ~1.9x that scope: **~$95-190K, ~19-57 channels** -- not the unscaled $50-100K/10-30 channels. For an explicit cell/channel budget from first principles instead of Replogle's leaner per-pert average, use `cells_needed = target_genes x cells_per_perturbation (500-1,000) / cells_per_channel (5,000-10,000)`; this DE-power-first route will exceed Replogle's own screened-scale numbers.

```python
# Replogle-style genome-wide design
# Each cell -> 1 library element (low MOI)
# Each gene -> 1 dual-sgRNA CRISPRi element (2 distinct sgRNAs per element)
# Replogle 2022 retained >2.5M cells at a median >100 cells per perturbation
# Total: ~9,900 expressed genes x 1 element = ~9,900 elements

```
