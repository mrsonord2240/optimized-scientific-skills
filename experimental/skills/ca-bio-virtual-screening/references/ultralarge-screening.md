# Hierarchical and ultralarge-library screening (moved from SKILL.md)

## Virtual Screening Pipeline (Hierarchical)

**Goal:** Screen 10M-compound library down to top-1k candidates for follow-up.

**Approach:** Three-stage filter. The 1% and top-1000 selections below are repository starting heuristics; choose production cutoffs from target-relevant enrichment, diversity, and throughput measurements.

Pseudo-code skeleton (orchestrator). Each helper function delegates to a dedicated skill: drug-likeness filter to `chemoinformatics/admet-prediction`, single-ligand Vina/GNINA to `dock_single` in `scripts/dock_single.py`, PoseBusters QC to `chemoinformatics/pose-validation`.

```python
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# Stub helpers to be implemented per project; see the cross-referenced skills.
def drug_like_filter(df):
    raise NotImplementedError('Implement via chemoinformatics/admet-prediction (Lipinski+Veber+PAINS)')
def vina_dock(smi, receptor_pdbqt, center, box):
    raise NotImplementedError('Wrap dock_single() from scripts/dock_single.py; return best affinity')
def gnina_rescore(smi, receptor_pdbqt, center, box):
    raise NotImplementedError('Wrap gnina --cnn_scoring rescore subprocess call')
def pose_validate(df):
    raise NotImplementedError('Implement via chemoinformatics/pose-validation (PoseBusters)')

def vs_pipeline(library_smi, receptor_pdbqt, center, box, output_dir, n_workers=16):
    df = pd.read_csv(library_smi)
    df_stage1 = drug_like_filter(df)

    worker = partial(vina_dock, receptor_pdbqt=receptor_pdbqt,
                     center=center, box=box)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        affinities = list(ex.map(worker, df_stage1['smiles']))
    df_stage1['vina_affinity'] = affinities
    df_stage2 = df_stage1.nsmallest(int(len(df_stage1) * 0.01), 'vina_affinity')

    df_stage2['gnina_affinity'] = df_stage2['smiles'].apply(
        lambda smi: gnina_rescore(smi, receptor_pdbqt, center, box))
    df_stage3 = df_stage2.nsmallest(1000, 'gnina_affinity')

    return pose_validate(df_stage3)
```

For very large libraries, use a restartable scheduler-backed workflow and measure throughput on a representative tranche. Record hardware, software version, box dimensions, ligand flexibility, and failure rate with every throughput estimate.

## Ultralarge Library Screening (ZINC22, Enamine REAL)

| Library | Scope | Typical access | Verification requirement |
|---------|-------|----------------|--------------------------|
| ZINC22 | Purchasable and make-on-demand compounds | Tranche/download interfaces | Record the tranche query and retrieval date |
| Enamine REAL | Make-on-demand compounds | Provider files or search interface | Record product-space release and retrieval date |
| Enamine HTS | Screening collection | Provider files | Confirm current stock/version with the provider |
| Mcule | Aggregated purchasable compounds | Provider search/export | Record filters and retrieval date |
| ChEMBL | Curated compounds and bioactivities | Versioned database release | Record ChEMBL release and extraction query |

Library sizes and availability change frequently. Obtain counts from the provider or versioned database at execution time rather than copying a static total into a workflow.

For ultralarge VS, the following percentages and thresholds are repository starting heuristics that must be calibrated for the target and library:
1. Apply a documented property/alert policy while retaining flagged and rejected counts
2. If known actives exist, test a permissive 2D-similarity prefilter such as ECFP4 Tanimoto >=0.4 and measure active/chemotype retention
3. Vina dock the filtered subset
4. Rescore top 1% with GNINA
5. Rescore top 0.1% with MM/GBSA or FEP

Lyu et al. (2019) screened 170 million make-on-demand compounds against AmpC and the D4 dopamine receptor. Of 549 D4 candidates synthesized and tested, 81 were new active chemotypes and 30 had submicromolar activity.
