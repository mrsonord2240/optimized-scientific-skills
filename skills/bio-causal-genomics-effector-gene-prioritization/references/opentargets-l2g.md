Moved out of `SKILL.md`: read when querying Open Targets L2G scores (legacy vs Platform schema, GWAS-type credible sets).

## Open Targets L2G via GraphQL

**Goal:** Query pre-computed L2G scores for a study and locus without re-running the integrative pipeline.

**Approach:** Query the Open Targets GraphQL endpoint with a study ID and lead variant; parse per-gene scores.

### Open Targets Platform vs Genetics Portal (2024 Consolidation)

In 2024 Open Targets Genetics was merged into the Open Targets Platform GraphQL API at `api.platform.opentargets.org/api/v4/graphql`. The legacy Genetics endpoint (`api.genetics.opentargets.org/graphql`) still responds but is deprecated; new pipelines should target the Platform API. The schema also changed: the Platform exposes `credibleSet(studyLocusId)` with `l2GPredictions` (target, score, SHAP per-feature explainability), whereas the legacy schema exposed `studyLocus2GeneTable` with `yProbaModel` and per-component sub-scores.

Legacy (Genetics, deprecated):

```graphql
query L2G_legacy($studyId: String!, $variantId: String!) {
  studyLocus2GeneTable(studyId: $studyId, variantId: $variantId) {
    rows {
      gene { symbol }
      yProbaModel
      yProbaDistance
      yProbaMolecularQTL
      hasColoc
    }
  }
}
```

Modern (Platform, recommended):

```graphql
query L2G_modern($studyId: String!) {
  credibleSet(studyLocusId: $studyId) {
    l2GPredictions {
      rows {
        target { approvedSymbol }
        score
        features { name value shapValue }
        shapBaseValue
      }
    }
  }
}
```

```python
import requests
import pandas as pd

resp = requests.post('https://api.platform.opentargets.org/api/v4/graphql',
                     json={'query': '...modern L2G query...',
                           'variables': {'studyId': 'GCST006464_locus_42'}})
preds = resp.json()['data']['credibleSet']['l2GPredictions']['rows']
l2g_df = pd.json_normalize(preds).sort_values('score', ascending=False)
```

The headline `score` (Platform) corresponds to `yProbaModel` (legacy). Platform `features[].shapValue` values replace the legacy `yProba*` sub-scores and explain what drove the prediction. Genes whose SHAP is dominated by the distance feature but minimal on QTL or chromatin features are distance-only candidates; trust the integrated `score` as primary.

**L2G is populated for GWAS-type credible sets, not molecular-QTL ones.** Querying a gene's own `credibleSets` (e.g. via `target(ensemblId)`) mostly returns eqtl/pqtl/sqtl-type loci with `l2GPredictions.count == 0` -- L2G is computed against GWAS-trait credible sets, at the trait's own lead locus. Use the top-level `credibleSets(studyTypes: [gwas], ...)` query to find a GWAS-type `studyLocusId` first (verified live, 2026-09-18). `examples/opentargets_l2g_query.py` is a runnable, dependency-free (stdlib `urllib`/`json` only) template implementing exactly this two-step lookup and the modern query above; run it directly against the live, unauthenticated API.
