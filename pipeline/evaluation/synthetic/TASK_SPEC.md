# Contract for building one synthetic eval task

You are building ONE task under `F:\Bio-Evals\synthetic\<task_dir>\`. Read first, in this order:
1. `F:\Bio-Evals\synthetic\README.md` (what these are and their limits)
2. `F:\Bio-Evals\synthetic\fine_mapping_effector\` (generate.py, solvers.py, eval.json) and `F:\Bio-Evals\synthetic\cis_mr_bcaa\` as exemplars of the shape and rigour expected
3. `F:\Bio-Evals\synthetic\common\grading.py` (how answers are graded, using LatchBio's own graders from `F:\Bio-Evals\latch-eval-tools`)
4. The graders themselves: `F:\Bio-Evals\latch-eval-tools\src\latch_eval_tools\graders\` (numeric.py, label_set.py, marker_gene.py, multiple_choice.py, composite.py). Know exactly what each accepts.

## Purpose

An LLM agent will later attempt the task in a network-less Linux sandbox (Python 3 + numpy/pandas/scipy/statsmodels/sklearn, R with many CRAN/Bioconductor packages, plink, bcftools, samtools). We are measuring whether reference Skills (documents of domain guidance) help it. So the task must have a **single correct, checkable answer** that a competent analyst reaches with the right method and that a common shortcut misses.

## Non-negotiable rules

- **Independence.** Design the trap from the published pitfall named in your brief and from general domain practice. Do NOT read any Skill (`F:\optimized-scientific-skills\skills\*`) or upstream bioSkills text to design the task; that would make the test circular. You may read Skills AFTER the design is frozen, only to fill `relevant_skills`.
- **Ground truth by construction.** Generate data from a known model. Keep the truth in `truth.json` / the eval graders, never inside `data/`.
- **`data/` holds only what a real requester would hand over**, in formats real tools read. No column, file or value names that hint at the answer or the trap (no `is_true`, `decoy`, `confounded`, etc.). Total under ~15 MB.
- **The task text is what the agent sees.** Written like a real request from a colleague. State every convention needed to make the answer unambiguous (thresholds, units, scale, what to exclude), but do NOT name the pitfall, the right method, or hint at what to avoid. End with the exact `<EVAL_ANSWER>` schema. Refer to files as `data/<name>`.
- **Deterministic.** Fixed seed; `python generate.py` reproduces `data/` and `eval.json` byte-for-byte.
- **No LLM or API calls, no internet.** Do not run `F:\Bio-Evals\runner\run.py`; I will run the models. Do not edit anything outside your own task dir. Do not edit `common/` or `runner/`.

## Files you must produce in `synthetic/<task_dir>/`

- `generate.py`: writes `data/`, `truth.json`, and `eval.json`.
- `eval.json`:
  ```json
  {"id": "<task_dir>", "replicates": null, "synthetic": true,
   "task": "<agent-visible text>", "graders": [<LatchBio grader specs>],
   "data_dir": "data", "relevant_skills": ["bio-..."], "design_source": "<one line: the pitfall and its source>"}
  ```
  `graders` use LatchBio grader types exactly as in the exemplars. If truth comes from the reference solver, have `solvers.py --write-eval` freeze it, as `cis_mr_bcaa` does.
- `solvers.py`: importable without side effects beyond loading data. Exposes `correct()` returning the answer dict, plus one `trap_*()` per named failure mode, each returning the answer a person taking that shortcut would give. Running `python solvers.py` prints a table (solver, answer, per-grader P/f, overall) and ends with `DISCRIMINATES (as designed)` (exit 0) or `!! DOES NOT DISCRIMINATE` (exit 1). Expected: `correct` passes; every named trap fails. Any trap that passes must be labelled `BLIND SPOT` and recorded, not hidden.
- `NOTES.md`: not agent-visible. Design, the pitfall and its source, what the truth is and how it was built, tolerance rationale, known blind spots, ways the task could be gamed or is ambiguous.

## Quality bar

- The reference `correct()` must reach the answer using ordinary, defensible analysis. Verify it recovers the constructed truth within your tolerances. Report by how much.
- Tolerances must be tight enough that the shortcut fails and wide enough that legitimate method variants pass (e.g. DESeq2 vs edgeR vs limma-voom; different LD tools). Say how you chose them.
- Check the shortcut fails by margin, not by a hair, and that the correct answer passes by margin.
- Think about how a mid-tier LLM would actually go about it and where it could go wrong for reasons unrelated to the skill you are testing (file-format friction, ambiguous wording). Remove that friction.
- Keep local runtime for `generate.py` and `solvers.py` under about 3 minutes. Local Python 3.14 has numpy, scipy, pandas, scikit-learn; `pip install` anything else you need locally (statsmodels, rdkit...). If you need bcftools/R, run them in the sandbox image: `docker run --rm --network none -v "<abs host dir>:/w" bioeval-env <cmd>` from a Python `subprocess` call (do not launch docker from Git Bash: it rewrites paths).
- `relevant_skills`: the Skill ids (directory names under `F:\optimized-scientific-skills\skills\`) a competent agent would plausibly open for this task. Used only to report uptake.

## Report back (under 200 words)

What you built; the correct solver's result vs constructed truth (numbers); each trap's outcome; blind spots; ambiguity risks; anything you had to change from the brief and why.
