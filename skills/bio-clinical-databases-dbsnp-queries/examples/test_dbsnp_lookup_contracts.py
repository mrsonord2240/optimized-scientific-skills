"""Regression checks for dbsnp_lookup.py's merge and batch contracts.

Usage: python test_dbsnp_lookup_contracts.py
Runs without network access; mocked services verify that a merged rsID is
annotated through its canonical ID, duplicate queries are removed, and the
not-found result has the documented status contract.
"""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("dbsnp_lookup.py")
SPEC = importlib.util.spec_from_file_location("dbsnp_lookup", MODULE_PATH)
dbsnp_lookup = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(dbsnp_lookup)


class _FakeVariantInfo:
    queries = None

    def getvariants(self, queries, fields):
        self.__class__.queries = list(queries)
        return [{
            "query": "rs429358",
            "_id": "chr19:g.45411941T>C",
            "gnomad_exome": {"af": {"af": 0.138498}},
            "gnomad_genome": {"af": {"af": 0.164436}},
            "clinvar": {"rcv": {"clinical_significance": "pathogenic"}},
        }]


class DbSnpLookupContractTests(unittest.TestCase):
    def test_not_found_has_status_contract(self):
        with patch.object(dbsnp_lookup, "refsnp", return_value=None):
            result = dbsnp_lookup.resolve_merge_chain("rs999999999999")
        self.assertEqual(result, {
            "status": "not_found",
            "final_rsid": "999999999999",
            "chain": ["999999999999"],
        })

    def test_batch_uses_canonical_query_once(self):
        def resolve(rsid):
            if rsid == "rs630496":
                return {"status": "resolved", "final_rsid": "429358", "chain": ["630496", "429358"]}
            return {"status": "resolved", "final_rsid": "429358", "chain": ["429358"]}

        alleles = [{"ref": "T", "alt": "C"}]
        with patch.object(dbsnp_lookup.myvariant, "MyVariantInfo", _FakeVariantInfo), \
             patch.object(dbsnp_lookup, "resolve_merge_chain", side_effect=resolve), \
             patch.object(dbsnp_lookup, "refsnp", return_value={}), \
             patch.object(dbsnp_lookup, "alleles_grch38", return_value=alleles):
            frame = dbsnp_lookup.batch_normalize_rsids(["rs630496", "rs429358", "rs429358"])

        self.assertEqual(_FakeVariantInfo.queries, ["rs429358"])
        self.assertEqual(frame["input_rsid"].tolist(), ["rs630496", "rs429358"])
        self.assertEqual(frame.loc[0, "canonical_rsid"], "429358")
        self.assertEqual(frame.loc[0, "myvariant_ids"], "chr19:g.45411941T>C")
        self.assertEqual(frame.loc[1, "gnomad_v2_exome_af"], "0.138498")


if __name__ == "__main__":
    unittest.main(verbosity=2)
