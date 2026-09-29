---
name: bio-clinical-databases-acmg-classification
description: Builds a non-diagnostic ACMG/AMP or AMP/ASCO/CAP evidence ledger using current ClinGen SVI methods, calibrated PP3/BP4 and PS3/BS3 thresholds, PVS1 prerequisite gates, and explicit somatic tiers. Use for training or qualified-review support, never stand-alone patient diagnosis or treatment decisions.
tool_type: python
primary_tool: requests
license: MIT
author: GPTomics
---

# ACMG/AMP Variant Classification

Apply the germline ACMG/AMP framework or the distinct AMP/ASCO/CAP somatic framework. Preserve an evidence trail; do not mix the two classification systems.

## Check compatibility and authority

The runnable example targets requests 2.31+ and the current public coordinate-based GeneBe endpoint. AutoPVS1 and InterVar are documented external tools, not shipped integrations.

Before using a runnable surface:

- inspect Python packages with `pip show <package>` and `help(module.function)`;
- inspect a CLI with `<tool> --version`;
- adapt examples to the installed API after `ImportError`, `AttributeError`, or `TypeError`, rather than retrying unchanged;
- query the current CSpec REST service and record the exact specification version and access date before applying generic rules.

Use this authority order for germline classification:

1. gene/disease-specific VCEP CSpec;
2. ClinGen SVI specifications;
3. generic ACMG/AMP 2015 as fallback.

Record the rule set, version or access date, transcript, disease mechanism, inheritance model, and evidence source used.

## Enforce the practice boundary

This Skill is an educational evidence-ledger aid, not a medical device or autonomous clinical classifier. Do not use its output for patient diagnosis, treatment selection, reproductive decisions, or stand-alone clinical reporting. A qualified variant scientist or clinical geneticist must review the disease association, transcript relevance, mechanism, VCEP specification, assay validity, evidence provenance, conflicts, and uncertainty before any clinical use. Do not place protected health information in public APIs or examples.

## Classify a germline variant

1. Normalize the variant and evaluate it on the clinically relevant transcript, preferring MANE Select where applicable.
2. Confirm the gene-disease relationship and mechanism. Do not apply PVS1 unless loss of function is an established mechanism for that disease.
3. Select the variant-type route:
   - predicted loss of function: apply the Abou Tayoun PVS1 decision tree and any VCEP-specific tree;
   - missense: use one calibrated predictor, not stacked correlated predictors;
   - splice or synonymous: evaluate SpliceAI and use transcript consequence or RNA evidence before escalating strength;
   - population evidence: compare gnomAD `grpmax_faf95` with BA1 and gene-specific Whiffin BS1 thresholds;
   - functional evidence: calibrate PS3/BS3 with the Brnich OddsPath framework;
   - segregation or in-trans evidence: apply the relevant ClinGen scoring guidance.
4. Apply current-rule integrity checks before summation: reject retired PP5/BP6, retain at most one strength or alias from each evidence family, surface opposed PP3/BP4 as a conflict rather than cancelling them, and do not count PP3 or PM4 in addition to PVS1 when they represent the same mechanism.
5. Convert each retained criterion to Tavtigian points, sum once, and classify:
   - Pathogenic: at least 10;
   - Likely Pathogenic: 6 to 9;
   - VUS: 0 to 5;
   - Likely Benign: -1 to -6;
   - Benign: -7 or less, or BA1 standalone when valid.
6. Report the criteria, strengths, evidence, overrides, exclusions, point total, provisional category, unresolved conflicts, and the required qualified-review/non-diagnostic boundary.

Read [the framework reference](references/acmg-framework.md) for calibrated thresholds, PVS1, splicing, population evidence, functional evidence, reconciliation, failure modes, and citations.

## Classify a somatic cancer variant

Use AMP/ASCO/CAP 2017 Tier I-IV, not germline P/LP/VUS/LB/B. Supply an explicit evidence category, same-tumor status, clinical significance, knowledgebase evidence record, and access date. Do not substitute a loose numeric OncoKB level. Tier III is uncertain significance; Tier IV is likely benign/benign. Cross-check current knowledgebases and retain the source record. See the somatic framework and cautions in [the framework reference](references/acmg-framework.md).

## Use the reusable example

Run [`scripts/acmg_classify.py`](scripts/acmg_classify.py) for the source example's scoring helpers and demonstration:

```bash
python scripts/acmg_classify.py
```

Its only third-party dependency is `requests`. The default demonstration is local, uses one computational predictor, prints unresolved context, and marks its category non-diagnostic. The Python file includes current coordinate-based GeneBe and versioned-gene CSpec helpers; use only public variant examples and never send patient data. PVS1 returns a review-required stop state when splice consequence or rescue-transcript review is missing. Every germline and somatic result retains the qualified-review boundary.

## Reference routes

| File | Read when |
|---|---|
| [`references/acmg-framework.md`](references/acmg-framework.md) | applying thresholds, decision trees, evidence strengths, failure modes, reconciliation, or reviewing citations |
| [`references/provenance.md`](references/provenance.md) | checking origin, attribution, or normalization scope |
| [`LICENSE`](LICENSE) | reviewing the retained MIT license terms |

## Related skills

- `bio-clinical-databases-variant-prioritization`: upstream rare-disease filtering.
- `bio-clinical-databases-clinvar-lookup`: ClinVar evidence aggregation.
- `bio-clinical-databases-gnomad-frequencies`: FAF95 support for BS1/BA1.
- `bio-clinical-databases-myvariant-queries`: aggregated annotations.
- `bio-variant-calling-clinical-interpretation`: clinical reporting workflow.
