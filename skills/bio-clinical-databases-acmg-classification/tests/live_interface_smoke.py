"""Bounded public-interface smoke using documentation examples only.

Never pass patient or private variant data to this script.
"""
import importlib.util
import pathlib


script = pathlib.Path(__file__).parents[1] / "scripts" / "acmg_classify.py"
spec = importlib.util.spec_from_file_location("acmg_classify", script)
acmg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(acmg)

genebe = acmg.genebe_api("17", 28364411, "C", "A")
variants = genebe.get("variants", [])
assert isinstance(variants, list) and variants, "GeneBe returned no variants"
print(f"GeneBe coordinate contract: PASS ({len(variants)} variant record(s))")

cspec = acmg.cspec_gene_versions("GATM")
versions = cspec.get("data", [])
assert isinstance(versions, list) and versions, "CSpec returned no version records"
assert all("/version/" in row.get("@id", "") for row in versions)
print(f"CSpec versioned-gene contract: PASS ({len(versions)} version record(s))")
