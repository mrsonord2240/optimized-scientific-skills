# Evaluate annotation evidence

```bash
python scripts/match_library.py references.mgf queries.mgf --tolerance 0.005 > matches.txt
```

This search supplies one evidence source; see [library inputs and limitations](library.md). Review exact mass, adduct, isotope pattern, spectrum, retention time and standard evidence as available, including unresolved isomers. A matching authentic standard under the same method provides stronger confirmation than a library score alone. Formula, class and candidate-set annotations are not equivalent to a confirmed structure. Molecular-network propagation is an analogue hypothesis; ambiguous candidates should not each count as independent confirmed compounds in enrichment. No automated comprehensive confidence adjudicator is supplied.

Done when each feature retains the evidence, alternatives and defensible confidence level with uncertainty.
