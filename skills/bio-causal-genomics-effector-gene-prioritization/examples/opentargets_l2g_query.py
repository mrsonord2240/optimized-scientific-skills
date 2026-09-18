#!/usr/bin/env python3
"""
Reference: Open Targets Platform GraphQL API (api.platform.opentargets.org/api/v4/graphql),
public, unauthenticated. Verified live 2026-09-18. Stdlib only (urllib + json) -- no
`requests`/`pandas` dependency required to run this.

Query pre-computed Open Targets L2G (locus-to-gene) scores for a GWAS credible set,
using the "modern" Platform query shape documented in this Skill's SKILL.md
(credibleSet -> l2GPredictions -> rows { target, score, features, shapBaseValue }).

Usage:
    python opentargets_l2g_query.py                      # demo: a real GWAS study locus
    python opentargets_l2g_query.py <studyLocusId>        # look up a specific study locus

Note: L2G is computed for GWAS-trait credible sets, not molecular-QTL (eqtl/pqtl/sqtl)
credible sets -- querying a gene's QTL-type credible sets (e.g. via `target(ensemblId)
{ credibleSets }`) will return l2GPredictions.count == 0 in the current Platform release.
Use `credibleSets(studyTypes: [gwas])` (as this script does) or a known GWAS studyLocusId.
"""
import json
import sys
import urllib.request

GQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"

FIND_GWAS_LOCUS_QUERY = """
query FindGwasLocus {
  credibleSets(studyTypes: [gwas], page: {index: 0, size: 5}) {
    count
    rows { studyLocusId studyId studyType }
  }
}
"""

L2G_QUERY = """
query L2G_modern($studyLocusId: String!) {
  credibleSet(studyLocusId: $studyLocusId) {
    studyLocusId
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
"""


def gql(query, variables=None):
    payload = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
    req = urllib.request.Request(
        GQL_URL, data=payload, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=30) as resp:
        body = json.load(resp)
    if "errors" in body:
        raise RuntimeError(f"GraphQL error: {body['errors']}")
    return body["data"]


def find_a_gwas_study_locus():
    """Demo helper: any real GWAS-type credible set currently in the Platform release."""
    data = gql(FIND_GWAS_LOCUS_QUERY)
    rows = data["credibleSets"]["rows"]
    if not rows:
        raise RuntimeError("No GWAS credible sets returned -- API schema may have changed.")
    return rows[0]["studyLocusId"]


def fetch_l2g(study_locus_id):
    data = gql(L2G_QUERY, {"studyLocusId": study_locus_id})
    credible_set = data.get("credibleSet")
    if credible_set is None:
        raise RuntimeError(f"studyLocusId {study_locus_id} not found.")
    return credible_set["l2GPredictions"]["rows"]


def main():
    study_locus_id = sys.argv[1] if len(sys.argv) > 1 else find_a_gwas_study_locus()
    print(f"studyLocusId: {study_locus_id}")

    rows = fetch_l2g(study_locus_id)
    if not rows:
        print("l2GPredictions: none for this credible set "
              "(common for molecular-QTL loci; try a GWAS-type studyLocusId).")
        return

    rows_sorted = sorted(rows, key=lambda r: r["score"], reverse=True)
    print(f"{'gene':<15}{'L2G score':>12}")
    for r in rows_sorted:
        print(f"{r['target']['approvedSymbol']:<15}{r['score']:>12.4f}")

    top = rows_sorted[0]
    distance_shap = sum(
        f["shapValue"] for f in top["features"] if f["name"].lower().startswith("distance")
    )
    qtl_shap = sum(
        f["shapValue"] for f in top["features"] if "qtl" in f["name"].lower()
    )
    print(f"\nTop gene {top['target']['approvedSymbol']} (score={top['score']:.4f}): "
          f"summed distance-feature SHAP={distance_shap:.4f}, summed QTL-coloc SHAP={qtl_shap:.4f}")
    if distance_shap > 0 and qtl_shap <= 0:
        print("-> distance-dominated call: treat as a weaker, distance-only candidate "
              "(see SKILL.md 'Nearest-gene assumption fails').")


if __name__ == "__main__":
    main()
