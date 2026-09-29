---
name: bio-clip-seq-ago-clip-mirna-targets
description: Identify direct miRNA-target interactions from AGO HITS-CLIP, AGO-CLEAR-CLIP (chimeric reads), HEAP (Halo-Ago2 mouse), chimeric eCLIP / miR-eCLIP (deep miRNA-target profiling), or CLASH using chimeric-read processing pipelines, seed-pairing analysis, and 3' auxiliary pairing rules. Use when distinguishing direct miRNA targets from indirect, integrating CLIP-derived target maps with TargetScan / miRDB / DIANA predictions, applying canonical 7mer-8mer seed matching with 3' UTR context, or recovering miRNA-mRNA chimeras at scale.
tool_type: mixed
primary_tool: chimeric-eCLIP
license: MIT
author: GPTomics
---

# AGO-CLIP and miRNA Target Identification

Use AGO CLIP-seq variants to map Argonaute-bound RNA and determine which
miRNA pairs with each target. Distinguish direct pairing evidence from
computational inference in every result:

- Chimeric methods such as CLEAR-CLIP, CLASH, and chimeric eCLIP / miR-eCLIP
  ligate the miRNA to its target and can identify the pair directly.
- Standard AGO HITS-CLIP, PAR-CLIP, and eCLIP identify AGO-bound sites but do
  not preserve miRNA identity; pair them with seed matching and matched miRNA
  expression data.
- HEAP profiles HaloTag-Ago2 in vivo in a transgenic mouse model and does not
  transfer directly to human samples.

Read [method selection and pairing rules](references/method-selection-and-pairing.md)
before choosing a workflow. Use [workflows and commands](references/workflows-and-commands.md)
for the source method patterns and [failure modes](references/failure-modes.md)
when results are sparse, ambiguous, or discordant.

## Verify the environment

The executable route is pinned to public Hyb commit
`028ab6371ce793ca5e86f475fce1f2cc6ad3c677`. Hyb is a Make-based program
without a semantic release identity; its `--version` reports GNU Make. The
Yeo chimeric-eCLIP workflow is a separate CWL project. HEAP is an experimental
method with public data and CLIPanalyze, not a standalone HEAP CLI. There is no
public `pyHyb 0.4+` distribution; `hybkit 0.3.6` is a distinct Hyb-format
library and is not an automatic replacement. TargetScanHuman 8, miRDB 6,
DIANA microT-CDS, samtools, and bedtools are separate prediction, service, or
command surfaces whose release and access status must be recorded.

Before using a code pattern, inspect the installed interface:

- Python: run `pip show <package>` and inspect the called function signature.
- CLI: run `<tool> --version` and `<tool> --help`.

If an installed interface differs, adapt the invocation rather than retrying
unchanged. Tooling workers must verify the exact Hyb interface and advertised
output schema before treating the supplied command patterns as executable.

## Choose the evidence route

| Goal | Route | Interpretation |
|---|---|---|
| Direct global miRNA-target pairs | Chimeric eCLIP, CLEAR-CLIP, or CLASH | Chimera is direct pairing evidence |
| Deep targets for one or a few miRNAs | miR-eCLIP with probe or PCR enrichment | Direct, enrichment-specific evidence |
| All AGO-bound sites | Standard AGO eCLIP / HITS-CLIP | Binding site only; miRNA identity remains inferred |
| In vivo mouse tissue | HEAP | Mouse-model evidence |
| Cost-conscious discovery | AGO HITS-CLIP plus TargetScan | Indirect computational assignment |
| Non-canonical or 3'-compensatory pairing | Chimeric method | Avoid a seed-only exclusion rule |

Do not report a per-miRNA direct target list from standard AGO-CLIP alone.
When the library is not chimeric, label assignments as predictions and filter
candidate miRNAs using matched small-RNA-seq expression.

## Run a chimeric analysis

1. Confirm the library method and whether miRNA-target ligation or probe
   enrichment occurred during library preparation.
2. Select preprocessing from the declared library layout. Total chimeric-eCLIP
   and targeted miR-eCLIP do not share a generic paired-UMI contract; follow the
   pinned protocol route in [workflows and commands](references/workflows-and-commands.md).
3. Prepare and install a named Hyb database under `HYB_HOME/data/db`, retaining
   the exact mature-miRNA and transcript reference versions and index manifest.
4. Identify miRNA-mRNA chimeras, retaining strand and alignment-quality
   information and separating contamination or non-mRNA partners.
5. Filter computational assignments using matched small-RNA-seq. Choose and
   justify a context-specific expression threshold, retain the values and
   units, and do not promote an inherited example threshold to a universal
   rule.
6. Aggregate stable read assignments while retaining site-level evidence.
   Report counts as assay recovery support, not binding affinity; normalize or
   stratify by declared expression, library depth, ligation, and UMI policy.
7. Cross-reference TargetScan conserved 7mer-m8 / 8mer predictions, but retain
   direct non-canonical chimeras as direct observations.
8. Validate top biological claims with an orthogonal reporter assay such as a
   luciferase or GFP fusion carrying the target 3' UTR.

Run [`scripts/run_chimeric_eclip.sh`](scripts/run_chimeric_eclip.sh) only after
the named database is installed. It executes two clean single-thread Hyb
replicates with explicit `detect`, `id`, `db`, and `type=mim` values, validates
the generated 16-column files, and retains only identical assignments present
in every replicate. It fixes the locale, hash seeds, Perl key perturbation, and
thread settings needed to stabilize the pinned source's count-tie traversal.
It publishes `sites.tsv`, `targets.tsv`, `support.tsv`, `excluded.tsv`, raw
logs, and `manifest.json` atomically. Review per-read support and exclusions;
do not reinterpret an unstable assignment as a direct target. When matched
expression filtering is enabled, expression values and the threshold must be
finite numbers; malformed or non-finite values stop the run with their input
location, and finite values below threshold receive the
`below_expression_threshold` exclusion reason.

## Run a computational assignment

1. Call AGO-bound peaks using the appropriate HITS-CLIP, PAR-CLIP, or eCLIP
   workflow.
2. Scan peak sequence or intersect annotated 3' UTR sites with TargetScan,
   miRDB, or DIANA candidates. TargetScan 8 coordinates are transcript/UTR
   relative, not genomic BED. Convert them with a version-matched spliced UTR
   map through `scripts/targetscan_sites_to_bed12.py`, then use `bedtools
   intersect -split -s` only when both datasets share assembly and strand.
3. Prioritize canonical 7mer-m8, 7mer-A1, and 8mer sites. Report 6mer evidence
   separately because it is common and weak.
4. Restrict candidate miRNAs to those expressed in the matched cell type and
   record the expression source and threshold.
5. Report `prediction + AGO overlap` as high-confidence indirect evidence, not
   as a recovered miRNA-mRNA chimera.
6. Route suspected non-canonical interactions to chimeric evidence or a
   separately validated full-duplex prediction workflow.

## Reconcile evidence

- A chimera without a TargetScan site may reflect non-canonical or
  3'-supplementary pairing; do not discard direct evidence solely for that
  reason.
- A TargetScan prediction with an AGO peak but no chimera may be functional,
  while chimera capture remains stochastic.
- A TargetScan prediction without an AGO peak is unsupported in that assay,
  tissue, condition, depth, and peak-calling context. It is not thereby proven
  biologically false or non-functional.
- An AGO peak without a TargetScan match remains an AGO-binding observation;
  investigate non-canonical pairing and expressed miRNAs before assigning it.
- Compare HEAP and human cell-line eCLIP only with explicit species and model
  context.

## Report

Record the library method, organism and tissue or cell type, AGO protein,
reference versions, chimera tool and mode, filtering thresholds, matched miRNA
expression source, seed classes, site and target counts, strand handling,
controls, and validation status. Separate direct chimeras, AGO-overlapped
predictions, and prediction-only candidates in the output. Include the exact
Hyb commit/database, deterministic process controls, consensus replicate
count, per-read assignment support, excluded-read reason counts, input/output
hashes, expression values and source, targeted-library UMI declaration and
protocol source when applicable, assembly, annotation and TargetScan releases,
and coordinate transformation manifest.

## Resources

- [Method selection and pairing rules](references/method-selection-and-pairing.md)
- [Workflows and command patterns](references/workflows-and-commands.md)
- [Failure modes and diagnostics](references/failure-modes.md)
- [Provider provenance](references/provenance.md)
- [Reusable chimeric eCLIP wrapper](scripts/run_chimeric_eclip.sh)
- [Orientation-aware consensus parser](scripts/consensus_hyb.py)
- [Protocol-declared targeted UMI extractor](scripts/extract_targeted_umi.py)
- [TargetScan UTR-to-genome BED12 converter](scripts/targetscan_sites_to_bed12.py)
