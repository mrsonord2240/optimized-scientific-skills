# Evaluation method and retained results

## Scope and status

The evaluation asked whether reference Skills help an agent solve bounded
scientific tasks, and whether packaging changes affect skill use. It began as
an upstream/optimized/no-Skills comparison, then focused on a selective package
and shorter descriptions after a full-catalog pilot had very low uptake.

These are synthetic local experiments, **not benchmarks.bio scores**. The
original public benchmark data were not obtained. Two tasks recreate the
structure of public variantbench problems using fabricated data; eight further
tasks exercise scientific pitfalls. The early plans in
[`evaluation/history/`](evaluation/history/) are retained as historical plans,
not claims that every proposed experiment was completed.

The source snapshot preserves ten task definitions and 263 top-level
`result.json` records across thirteen run folders. Interrupted attempts remain
in the local evidence archive and are not counted as completed results.
Infrastructure failures remain visible and are excluded from model pass-rate
denominators. Current harness files are the final local snapshot; they are not
proof that every older run used identical code or environment bytes.

## Task and execution contract

- Generate a synthetic world with a fixed seed and a checkable answer. Save the
  task, answer key, generator, reference solver, shortcut solvers, and known
  blind spots. A passing reference solver establishes the grader wiring, not
  the task's difficulty for an agent.
- Give the agent only the `task` string and `data/`. Keep `graders`, truth,
  generators, solvers, and design notes outside the sandbox.
- Use a bash-only interaction loop. The agent emits one shell block per turn
  and ends with JSON inside `<EVAL_ANSWER>` tags. Grade on the host.
- Run the analysis in a network-disabled Docker container. The final runner
  copies inputs into the container and uses 2 CPUs and 8 GiB of memory.
- Change only the skill arm and catalog presentation within a comparison.
  `none` has no Skills; `upstream` uses original Skills; `optimized` uses the
  published shelf; `draft` contains uncertified experimental edits.
- For `--bundle domain`, deliver twelve Skills: relevant Skills, siblings from
  their upstream domains, and deterministic task-seeded filler. Full
  descriptions come from YAML frontmatter. `--desc short` extracts the first
  matching "Use ..." sentence and caps it at fifteen words.
- Record the answer, grader outcomes, finish reason, configuration, delivered
  Skills, observed skill references, cost, timing, and infrastructure status.
  Retain trajectories and events for inspection. The uptake metric observes
  skill paths in shell commands; it is a proxy for reading, not proof of use.

The final runner defaults to 60 steps, a 25-minute wall limit, a 300-second
command timeout, a 32,000-token reply cap, and a USD 1.50 cost stop. Use each
result's saved `config` and run label when comparing history; older pilots and
provider routes differed. Model strings below are historical recorded IDs.

## The completed selective-package comparison

Model: `deepseek/deepseek-v4.1-flash`. Seven tasks, five repeats per condition.
The saved arm manifest pins 183 source Skills to optimized commit
`2f381782596c6569fe5a8357556856512b7fbb6c` and upstream commit
`d91ed3d563019e649dc854c56ccd62551359488a`; only twelve were delivered per task.

| Condition | Passed | Observed any Skill | Infrastructure failures |
| --- | ---: | ---: | ---: |
| No Skills | 18/35 (51.4%) | 0/35 | 0 |
| Optimized, full descriptions | 19/35 (54.3%) | 17/35 (48.6%) | 0 |
| Optimized, short descriptions | 24/35 (68.6%) | 26/35 (74.3%) | 0 |

Each cell below is passes out of five:

| Task | No Skills | Full | Short |
| --- | ---: | ---: | ---: |
| `cis_mr_bcaa` | 4 | 3 | 3 |
| `fine_mapping_effector` | 2 | 1 | 5 |
| `i1_ml_leakage` | 3 | 2 | 2 |
| `i2_pseudobulk_de` | 1 | 1 | 2 |
| `i5_crispr_cn_bias` | 4 | 5 | 5 |
| `i6_proteomics_de` | 0 | 2 | 2 |
| `i8_microbiome_da` | 4 | 5 | 5 |

The baseline folder contains fifty results over ten tasks. The three extra
tasks (`i3_ora_background`, `i7_isotope_correction`, `i9_vcf_normalization`)
each passed 5/5 without Skills and are outside the seven-task comparison.
The summary script selects the common task set before counting the baseline;
comparing its unfiltered 33/50 to the two 35-run arms would be misleading.

This supports a descriptive observation: in this small selected task set,
short descriptions coincided with more observed skill references and more
passes. It does not establish a general causal effect. Four of the five
additional passes over full descriptions came from `fine_mapping_effector`.
There was no upstream arm at the twelve-Skill package size, so this experiment
does **not** establish that optimized Skills beat the originals. Nor does it
isolate the value of an uncapped "Use when" sentence from the fifteen-word cap.

## Other retained experiments

[`evaluation/summary.json`](evaluation/summary.json) lists every run folder,
model, arm, description mode, task, pass count, uptake, and infrastructure
failure count without pooling incompatible conditions.

- The 183-Skill pilot has six completed results across upstream and optimized
  arms. Neither arm recorded skill uptake. It is retained as a failed
  packaging experiment, not a powered upstream comparison.
- `runs_draft_short`, `runs_recut_short`, and `runs_recut_prot` test script and
  router changes after inspecting earlier failures. They are development
  evidence on reused tasks, not independent held-out evaluation. The later
  recuts passed pseudobulk 5/5, proteomics 5/5, and ML leakage 3/5.
- `runs_cli_opus` has 45 passes in 45 runs across three arms on three tasks.
  That ceiling provides no pass-rate discrimination between the arms.
- The two `runs_strong_*` folders contain thirteen infrastructure failures
  among 36 attempts. Keep those separate from task failures and from the CLI
  experiment. CLI cost fields are list-price estimates, not metered charges.
- `runs`, `runs_oldrule`, and `runs_base_v1_uncapped` retain older pilots and
  changed-rule/budget histories; do not combine them with the completed
  comparison as extra repeats.

## Reproduce or inspect

The compact checkout supports offline recounting with standard-library Python:

```text
python pipeline/evaluation/summarize.py --check
python pipeline/evaluation/summarize.py
```

The second command prints a deterministic JSON summary. It reads saved results
only and never invokes a model, Docker, a solver, or the historical runner.
Every imported result records its source path and original file SHA-256 in
[`evaluation/results.jsonl`](evaluation/results.jsonl).

To restore the historical working materials, use the local archive inventory
in [EVIDENCE.md](EVIDENCE.md). It includes task data, arm bytes (including draft
edits), run results, trajectories, events, and build logs. Credentials, Python
caches, and redundant per-run workspace copies are excluded. Restore into a
new directory so the original evidence stays intact.

For an independently authorized future rerun:

1. Restore task inputs and the exact arm snapshot, or regenerate inputs in a
   separate workspace and compare their hashes with the archive inventory.
   The existing `build_arms.py` is preserved for provenance: it has a hardcoded
   shelf path and replaces arm directories. Do not use it to reconstruct an
   old arm from today's `HEAD`.
2. Obtain `latchbio/latch-eval-tools` at
   `75cd1e8fc471e5622ec020a9e3bbf604ce88b454` beside `synthetic/` and `runner/`.
   Its grader code is not redistributed here; its checked-in license says
   "All rights reserved." The adapter loads the grader modules without the
   full Latch SDK. The task generators also require NumPy, SciPy, pandas, and
   task-specific packages documented in their source and notes.
3. Inspect/build the preserved Dockerfile, record the resulting image digest
   and installed versions, and verify required tools. The historical image
   was not pinned by digest; several installs deliberately tolerate failure.
   Rebuilding it today cannot promise the old environment. The final runner
   also requires `requests` and PyYAML on its host.
4. Check the reference and shortcut solvers against the exact inputs before
   spending on models. Keep keys outside the repository. Use explicit task,
   arm, model, repeat count, description mode, and fresh output directory.
   Example runner flags are in `runner/run.py`; defaults include all arms.
5. Preserve complete trajectories and failures; summarize matched tasks and
   repeats. Never treat post-failure development runs as held-out evidence or
   substitute a successful exit for scientifically correct output.

This closeout validated preservation, hashes, syntax, links, and recounts.
It did not rerun the scientific tasks or certify the historical environment.
