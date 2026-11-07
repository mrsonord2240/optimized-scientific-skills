# Readiness and records

Use the current `skill-auditor.zip` in the records repository as the canonical
rubric and report-schema source. Extract it only to a disposable run directory
when it is not already installed.

Unless the current rubric is stricter, `Production Ready` requires final score
at least 85, static score at least 80, execution average at least 85, Layer 1
at least 32, Layer 2 at least 48, assertion pass rate at least 90 percent, no
veto, no open P0, and representative inspected execution for every materially
distinct required workflow and output family. Exact-identity reusable evidence
may satisfy unchanged surfaces under the re-audit contract.

The Marketplace does not read audit scores: providers do not submit evaluation
reports, and its `intake:skill` gate checks structure, metadata, and license
evidence only. The audit is our evidence that the science and commands hold,
for the maintainer who decides inclusion. Treat the gate above as the finish
line; a higher score is not a goal.

Report recommendations use `P0`, `P1`, or `P2` only. The records indexer
rejects any other priority, so map a would-be P3 to P2 or drop it.

Use these states:

- `untouched`: no usable audit for the relevant source;
- `audited`: audit exists but final certification is incomplete;
- `candidate-ready`: independent risk-based final audit passes for exact
  working bytes with representative coverage of every materially distinct
  required workflow and output family, but they have not entered the
  run-closing product commit;
- `ready`: exact audited bytes are committed on the optimized shelf with
  reconciled provider metadata;
- `done`: Marketplace intake accepted the exact optimized commit and path.

For every audit identity, preserve schema-valid `report.json`, `viewer.md`,
saved run scripts and inputs, command/output evidence, source identity, and
explicit execution classifications. Do not claim evidence that exists only in
chat.

Publish the initial baseline and accepted final audit with the repository
tooling. Rejected intermediate runs may remain in durable local evidence unless
publication is needed to preserve a blocker, explain an identity transition, or
satisfy the records schema. Legacy flat audit roots remain supported:

```powershell
python tools/publish_audits.py --repo <provider-root> --skill <skill-id>
npm run audits:index
```

For the active modular run layout, preserve the strict report unchanged and
pass source/candidate provenance through `source-identity.json`, which must
record the candidate's absolute `path` (candidate binding reads it). Each run
needs `report.json`, `viewer.md`, and `source-identity.json` at the run root;
artifacts must sit under the run's `scripts/`. The publisher refuses to
overwrite a published version, so every audit run uses a new run directory
name. Publish the exact run explicitly and enumerate every saved script or
bounded input selected for the public record:

```powershell
python tools/publish_audits.py --repo <records-root> --skill <skill-id> `
  --run-dir <raw-run-root> --artifact <script-or-input> [...]
npm run audits:index
npm run audits:check
```

Candidate records use a deterministic content or Git-tree identity and retain
the immutable origin for attribution. They do not satisfy provider-promotion
gates until exact audited bytes are committed on the optimized shelf and the
final accepted audit is bound to that full provider commit.

Regenerate the canonical index, backlog, `audits/STATUS.md`, and
`audits/STATUS.html` from records. The two status outputs present the same
canonical known Skill inventory by category and show at least total known,
audited, untouched, ready, and out-of-scope counts plus totals. `Audited`
includes ready Skills. Derive:

- category, total known, ready, and out-of-scope states from
  `audits/CORPUS.json`;
- audited from the latest published audit records;
- untouched as in-scope known Skills with no published audit.

Never edit counts by hand. Run `npm run audits:index` after every published
audit.

Treat provider inventory refresh as a mandatory post-commit invariant. After
every commit, amendment, or replacement commit that adds, changes, or removes a
Skill or its provider metadata in `optimized-scientific-skills`, run these from
the records repository before continuing:

```powershell
npm run audits:inventory
npm run audits:check
```

`audits:inventory` is the repository alias for refreshing `audits/CORPUS.json`
from the canonical optimized provider and regenerating every audit view. Its
expanded local form is:

```powershell
npm run audits:index -- --refresh-inventory --provider F:\optimized-scientific-skills
```

A plain `npm run audits:index` only regenerates views from the saved corpus and
does not satisfy this post-commit requirement. Do not accept the product commit
as reconciled until inventory refresh and consistency check both pass. If either
command fails, resolve or report the concrete blocker before accepting the
transition.
