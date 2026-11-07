# Execution and evidence standard

Save the runnable command or script, input identity, environment and relevant
versions, output location, exit status, and assertions. Use a bounded real
public input when available. Label synthetic fixtures and their planted truth.

Compilation, parsing, imports, `--help`, file existence, nonzero file size, and
exit code zero are preflight facts, not final behavioral evidence. Do not mark
an advertised path executed unless its code actually ran.

Match assertions to the artifact:

- parse structured data and check schema, dimensions, counts, ranges, units,
  and meaningful values;
- check scientific invariants, cohorts, denominators, transformations, and
  planted or independently expected relationships;
- render or open figures, tables, and reports at intended size and inspect
  labels, clipping, legibility, and interpretability;
- verify interactive outputs expose the expected information;
- exercise meaningful failure guards where the documented contract needs them.

An empty, malformed, misleading, illegible, scientifically unexpected, or
unusable output fails even if the process exits successfully.

Classify each surface as `executed`, `failed`, `static-only`, `blocked`, or
`not-applicable`. A heavy optional surface (see the tooling worker Skill) is
`static-only` with the reason `heavy-optional`; it does not block readiness
when the Skill labels it as not executed. Documentation and source inspection can support a finding but
are not executed evidence. A blocked route retains its limitation and
after-action item; it does not receive an inflated score.
