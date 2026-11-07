# Column statistics

```bash
python scripts/statistics.py columns alignment.fasta columns.json
```

Conservation ignores gaps and defaults to minimum occupancy 0.5; below it the score is null and excluded from the mean. Entropy is in bits; information content is KL divergence against the selected background. The source chooses uniform DNA or Robinson protein background from the sequence; inspect the reported label and excluded alphabet before interpreting results. For JSD functional-site scores read `references/capra-singh-jsd.md`; for entropy, PSSM, Neff and MI-APC read `references/information-content-pssm.md`; for gaps read `references/gap-statistics.md`.

For a centered moving-average conservation profile use `examples/conservation_profile.py` and its `conservation_profile` function; the source window default is 10. JSD uses a different background than the original published implementation and omits its sequence weights; see the JSD reference before interpreting rankings.

Done when requested columns or summary are reported with the denominator, occupancy and background.
