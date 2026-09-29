"""Focused tests for current PAML protein and codon marginal rst rows."""

import argparse
import importlib.util
import math
import pathlib
import re
import sys
import tempfile
import unittest


SKILL_ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = SKILL_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from paml_rst import parse_marginal_rst
import ancestral_reconstruction
import codeml_asr


HEADER = """PAML test output
9..10 10..1
Nodes 9 to 10 are ancestral
(1) Marginal reconstruction of ancestral sequences
Prob of best state at each node, listed by site
   site   Freq   Data:
"""


def parse_text(rows):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".rst", encoding="utf-8") as handle:
        handle.write(HEADER + rows)
        handle.flush()
        return parse_marginal_rst(handle.name)


class PamlRstTests(unittest.TestCase):
    def test_protein_matrix_preserves_legacy_state_and_probability(self):
        sequences, probabilities = parse_text(
            "   1      1   A-G : A(0.900) G(0.700)\n"
            "   2      1   A-G : C(0.800) G(1.000)\n"
        )
        self.assertEqual(sequences, {"Node10": "AC", "Node9": "GG"})
        self.assertEqual(probabilities["Node10"][0], {
            "node": "Node10", "site": 1, "state": "A", "probability": 0.9,
        })

    def test_codon_matrix_uses_codon_probability_and_retains_aa_probability(self):
        sequences, probabilities = parse_text(
            "   1      1   ATG : ATG M 0.100 (M 0.348) GGG G 0.800 (G 0.850)\n"
            "   2      1   TAT : TAT Y 0.100 (A 0.246) GAT D 0.500 (D 0.700)\n"
        )
        self.assertEqual(sequences, {"Node10": "ATGTAT", "Node9": "GGGGAT"})
        self.assertEqual(probabilities["Node10"][0], {
            "node": "Node10", "site": 1, "state": "ATG", "probability": 0.1,
            "codon_translation": "M", "amino_acid_state": "M", "amino_acid_probability": 0.348,
        })
        self.assertEqual(probabilities["Node10"][1]["codon_translation"], "Y")
        self.assertEqual(probabilities["Node10"][1]["amino_acid_state"], "A")
        self.assertEqual(probabilities["Node10"][1]["amino_acid_probability"], 0.246)

    def test_no_data_and_malformed_rows_remain_explicit(self):
        with self.assertRaisesRegex(ValueError, "no marginal site rows"):
            parse_text("")
        with self.assertRaisesRegex(ValueError, "malformed codon node values"):
            parse_text("   1      1   ATG: ATG M 1.100 (M 1.100) GGG G 0.8 (G 0.8)\n")
        self.assertEqual(ancestral_reconstruction.assess_reconstruction_quality([])["status"], "no_data")

    def test_codeml_wrapper_exposes_codon_and_amino_acid_posteriors(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".rst", encoding="utf-8") as handle:
            handle.write(HEADER + "   1      1   ATG : ATG M 0.100 (M 0.348) GGG G 0.800 (G 0.850)\n")
            handle.flush()
            parsed = codeml_asr.parse_rst_posteriors(handle.name)
        self.assertEqual(parsed[10][0], {
            "site": 1, "state": "ATG", "prob": 0.1,
            "codon_translation": "M", "amino_acid_state": "M", "amino_acid_probability": 0.348,
        })


def parse_real_source_rows(path, expected_sites, codon=False):
    """Read PAML's printed best-state values as the source-side oracle."""
    row_pattern = re.compile(r"^\s*(\d+)\s+\d+\s+.*?:\s*(.*)$")
    value_pattern = (
        re.compile(r"([ACGTN?\-]{3})\s+([A-Za-z*?\-])\s+([0-9.]+)\s+"
                   r"\(([A-Za-z*?\-])\s+([0-9.]+)\)")
        if codon else re.compile(r"([A-Za-z*?\-])\(([0-9.]+)\)")
    )
    rows = {}
    in_rows = False
    for line in pathlib.Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if "site" in line.lower() and "freq" in line.lower() and "data:" in line.lower():
            in_rows = True
            continue
        if not in_rows:
            continue
        match = row_pattern.match(line)
        if not match:
            continue
        site, body = int(match.group(1)), match.group(2)
        matches = list(value_pattern.finditer(body))
        if not matches or "".join(item.group(0) for item in matches).replace(" ", "") == "":
            raise AssertionError(f"{path}: could not read source values for site {site}")
        values = []
        for item in matches:
            if codon:
                codon_state, translation, probability, aa_state, aa_probability = item.groups()
                values.append((codon_state.upper(), float(probability), translation.upper(),
                               aa_state.upper(), float(aa_probability)))
            else:
                state, probability = item.groups()
                values.append((state.upper(), float(probability)))
        rows[site] = values
        if len(rows) == expected_sites:
            break
    if list(rows) != list(range(1, expected_sites + 1)):
        raise AssertionError(f"{path}: source rows do not cover sites 1..{expected_sites}")
    return rows


def assert_real_rst(path, expected_sites=110, state_width=1, codon=False):
    ancestors, posteriors = parse_marginal_rst(path)
    if len(ancestors) != 7 or set(map(len, ancestors.values())) != {expected_sites * state_width}:
        raise AssertionError(f"{path}: expected 7 ancestors x {expected_sites} sites")
    if len(posteriors) != 7 or any(len(rows) != expected_sites for rows in posteriors.values()):
        raise AssertionError(f"{path}: posterior coverage does not match sequence coverage")
    for rows in posteriors.values():
        if [row["site"] for row in rows] != list(range(1, expected_sites + 1)):
            raise AssertionError(f"{path}: site indices are not 1..{expected_sites}")
        if any(not math.isfinite(row["probability"]) or not 0 <= row["probability"] <= 1 for row in rows):
            raise AssertionError(f"{path}: invalid best-state posterior")
    # Compare every parsed best state and probability with the values printed in
    # this exact successful PAML run. Posterior magnitudes can vary across valid
    # runs, so the fixture contract follows its own source output.
    raw_rows = parse_real_source_rows(path, expected_sites, codon=codon)
    node_order = list(posteriors)
    for site in range(1, expected_sites + 1):
        if len(raw_rows[site]) != len(node_order):
            raise AssertionError(f"{path}: source row {site} has unexpected node coverage")
        for node, source in zip(node_order, raw_rows[site]):
            parsed = posteriors[node][site - 1]
            if codon:
                actual = (parsed["state"], parsed["probability"], parsed["codon_translation"],
                          parsed["amino_acid_state"], parsed["amino_acid_probability"])
            else:
                actual = (parsed["state"], parsed["probability"])
            if actual != source:
                raise AssertionError(f"{path}: parsed posterior disagrees with PAML site {site}, {node}")
    return ancestors, posteriors


def run_real_fixture_checks(args):
    if args.codon_rst:
        ancestors, posteriors = assert_real_rst(args.codon_rst, state_width=3, codon=True)
        wrapped = codeml_asr.parse_rst_posteriors(args.codon_rst)
        assert len(wrapped) == 7 and len(wrapped[10]) == 110
        assert wrapped[10][1]["state"] == posteriors["Node10"][1]["state"]
        assert wrapped[10][1]["prob"] == posteriors["Node10"][1]["probability"]
        assert wrapped[10][1]["amino_acid_probability"] == posteriors["Node10"][1]["amino_acid_probability"]
        sample = posteriors["Node10"][:2]
        print(f"PASS real CODONML: {len(ancestors)} nodes x 110 sites; 770 codon and amino-acid posteriors match PAML source rows; Node10 sites1-2={[(x['state'], x['probability'], x['codon_translation'], x['amino_acid_state'], x['amino_acid_probability']) for x in sample]}")

    for label, path in (
        ("provider-protein", args.protein_provider_rst),
        ("extracted-protein", args.protein_extracted_rst),
    ):
        if not path:
            continue
        ancestors, posteriors = assert_real_rst(path)
        if label == "provider-protein":
            public_sequences = ancestral_reconstruction.parse_rst_ancestors(path)
            public_probabilities = ancestral_reconstruction.extract_site_probabilities(path)
            assert public_sequences == ancestors and len(public_probabilities) == 770
            assert public_probabilities[0]["state"] == posteriors["Node10"][0]["state"]
            assert public_probabilities[0]["probability"] == posteriors["Node10"][0]["probability"]
        else:
            wrapped = codeml_asr.parse_rst_posteriors(path)
            assert len(wrapped) == 7 and len(wrapped[10]) == 110
            assert wrapped[10][0]["state"] == posteriors["Node10"][0]["state"]
            assert wrapped[10][0]["prob"] == posteriors["Node10"][0]["probability"]
        sample = [(x["state"], x["probability"]) for x in posteriors["Node10"][:2]]
        print(f"PASS real {label}: {len(ancestors)} nodes x 110 sites; all parsed states/posteriors match PAML source rows; Node10 sites1-2={sample}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--codon-rst")
    parser.add_argument("--protein-provider-rst")
    parser.add_argument("--protein-extracted-rst")
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PamlRstTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    run_real_fixture_checks(args)


if __name__ == "__main__":
    main()
