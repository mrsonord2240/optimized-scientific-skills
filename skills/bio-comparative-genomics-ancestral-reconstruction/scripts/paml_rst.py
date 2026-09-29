"""Parsers for the matrix-oriented marginal reconstruction section in PAML rst."""

import re


_NUMBER = r"(?:0(?:\.\d+)?|1(?:\.0+)?)"
_PROTEIN_STATE_PROB = re.compile(rf"([A-Za-z*?\-])\(({_NUMBER})\)")
_CODON_STATE_PROB = re.compile(
    rf"([ACGTN?\-]{{3}})\s+([A-Za-z*?\-])\s+({_NUMBER})\s+"
    rf"\(([A-Za-z*?\-])\s+({_NUMBER})\)"
)


def _parse_node_values(text, site):
    """Parse one PAML node vector, retaining codon and amino-acid posteriors."""
    codon_matches = list(_CODON_STATE_PROB.finditer(text))
    if codon_matches:
        cursor = 0
        values = []
        for match in codon_matches:
            if text[cursor:match.start()].strip():
                raise ValueError(f"PAML rst site {site} has malformed codon node values")
            codon, amino_acid, probability, listed_amino_acid, amino_probability = match.groups()
            probability = float(probability)
            amino_probability = float(amino_probability)
            if probability > amino_probability + 0.002:
                raise ValueError(f"PAML rst site {site} codon posterior exceeds its amino-acid posterior")
            values.append((codon.upper(), probability, amino_acid.upper(),
                           listed_amino_acid.upper(), amino_probability))
            cursor = match.end()
        if text[cursor:].strip():
            raise ValueError(f"PAML rst site {site} has malformed codon node values")
        return "codon", values

    protein_matches = list(_PROTEIN_STATE_PROB.finditer(text))
    if not protein_matches:
        return None, []
    cursor = 0
    values = []
    for match in protein_matches:
        if text[cursor:match.start()].strip():
            raise ValueError(f"PAML rst site {site} has malformed protein node values")
        state, probability = match.groups()
        values.append((state.upper(), float(probability), None, None, None))
        cursor = match.end()
    if text[cursor:].strip():
        raise ValueError(f"PAML rst site {site} has malformed protein node values")
    return "protein", values


def parse_marginal_rst(rst_file, expected_nodes=None, expected_sites=None):
    """Return best-state sequences and probabilities keyed by PAML node number.

    CODEML's marginal section stores a site row followed by one best-state
    probability per reconstructed node. The node order is recoverable from the
    branch list and ancestral-node range in current PAML 4.10 output.
    """
    with open(rst_file, encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    ancestral = set()
    for line in lines:
        match = re.search(r"Nodes\s+(\d+)\s+to\s+(\d+)\s+are ancestral", line, re.I)
        if match:
            low, high = map(int, match.groups())
            ancestral.update(range(low, high + 1))
            break

    node_order = []
    marginal_start = next((i for i, line in enumerate(lines)
                           if "Prob of best state at each node" in line), len(lines))
    for line in lines[:marginal_start]:
        edge_pairs = re.findall(r"\b(\d+)\.\.\s*(\d+)\b", line)
        candidate = []
        for _, child in edge_pairs:
            child = int(child)
            if child in ancestral and child not in candidate:
                candidate.append(child)
        missing_root = ancestral.difference(candidate)
        if len(candidate) + len(missing_root) == len(ancestral) and len(missing_root) == 1:
            # PAML's edge list names descendants; the root is never a child.
            node_order = candidate + sorted(missing_root)
            break

    if not ancestral or not node_order:
        raise ValueError("PAML rst lacks the ancestral-node range or branch-order header")

    sequences = {node: [] for node in node_order}
    posteriors = {node: [] for node in node_order}
    site_ids = []
    row_format = None
    in_rows = False
    for line in lines:
        if "site" in line.lower() and "freq" in line.lower() and "data:" in line.lower():
            in_rows = True
            continue
        if not in_rows:
            continue
        match = re.match(r"\s*(\d+)\s+\d+\s+.*?:\s*(.*)$", line)
        if not match:
            if site_ids and ("Summary of changes" in line or "TREE #" in line):
                break
            continue
        site = int(match.group(1))
        current_format, values = _parse_node_values(match.group(2), site)
        if current_format is None:
            raise ValueError(f"PAML rst site {site} has malformed marginal node values")
        if row_format is None:
            row_format = current_format
        elif row_format != current_format:
            raise ValueError("PAML rst marginal rows mix protein and codon formats")
        if len(values) != len(node_order):
            raise ValueError(
                f"PAML rst site {site} has {len(values)} node values; expected {len(node_order)}"
            )
        site_ids.append(site)
        for node, (state, probability, codon_translation, amino_acid, amino_probability) in zip(node_order, values):
            probability = float(probability)
            if not 0.0 <= probability <= 1.0:
                raise ValueError(f"PAML rst site {site}, node {node} has invalid probability")
            sequences[node].append(state)
            record = {
                "node": f"Node{node}", "site": site,
                "state": state, "probability": probability,
            }
            if codon_translation is not None:
                record["codon_translation"] = codon_translation
                record["amino_acid_state"] = amino_acid
                record["amino_acid_probability"] = amino_probability
            posteriors[node].append(record)

    if not site_ids:
        raise ValueError("PAML rst contains no marginal site rows")
    if expected_nodes is not None and len(node_order) != expected_nodes:
        raise ValueError(f"PAML rst has {len(node_order)} ancestral nodes; expected {expected_nodes}")
    if expected_sites is not None and len(site_ids) != expected_sites:
        raise ValueError(f"PAML rst has {len(site_ids)} sites; expected {expected_sites}")
    if len(set(site_ids)) != len(site_ids) or site_ids != list(range(site_ids[0], site_ids[0] + len(site_ids))):
        raise ValueError("PAML rst marginal site indices are duplicated or non-contiguous")

    return ({f"Node{node}": "".join(sequences[node]) for node in node_order},
            {f"Node{node}": posteriors[node] for node in node_order})
