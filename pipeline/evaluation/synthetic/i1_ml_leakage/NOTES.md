# i1_ml_leakage (not agent-visible)

All data fabricated. `python generate.py` (~5 s) then `python solvers.py` (~15 s).

## Pitfall and source
(a) Feature selection done on all samples before cross-validation gives optimistic error on high-dimensional noise (Ambroise & McLachlan 2002, PNAS 99:6562). (b) Replicate samples from one patient split across folds leak patient identity (record-wise vs subject-wise CV). The agent is asked for the honest expected ROC AUC of "top-50 t-test -> z-score -> L2 logistic regression" on samples/patients not used for any fitting or selection. The task text never says leakage, nesting or grouping.

## Construction and truth
- **A**: 60 x 5000, iid N(0,1), balanced random labels. True AUC on new samples = 0.5 exactly, for any pipeline.
- **B**: 40 patients x 3 replicates x 2000 features, label per patient (20/20). 8 signal features (patient-level shift +-1.1, patient SD 0.5, replicate SD 1.0); 1200 fingerprint features (patient SD 1.2, replicate SD 0.4, ICC ~0.9, unrelated to label); 792 pure noise (SD 1). Columns and rows shuffled.
- **B truth = 0.7607** (graded as 0.76): the specified pipeline fitted on the released (4-dp rounded) data, scored on 4000 freshly drawn patients from the generating model. Independent of any CV code.
- Unconditional expectation over redrawn datasets is ~0.76 (scan of 20 worlds: mean 0.757, sd 0.028), so the realised dataset is typical.
- Leaky answers: selection-before-CV on A = 0.997; random-split or selection-outside on B = 1.00. Resubstitution = 1.00.

## Seeds (a constraint on the fabricated world)
A and B use independent RNG streams. `SEED_A=17`, `SEED_B=23` were chosen from scans of 1-30 / 1-24 so that the common honest CV variants (3/5/10-fold, LOO) land close to the truth on the released data. Honest CV on 60 noise samples is very noisy: across seeds, pooled 5-fold AUC on pure noise ranged 0.35-0.66 and LOO 0.13-0.75. Without this selection a legitimate analyst could land far from 0.5 on A by luck.

## Tolerances (numeric_tolerance, absolute, asymmetric; one grader per dataset)
- A: 0.5 -0.20/+0.20 -> [0.30, 0.70]. Wide because honest estimates on n=60 noise have SD ~0.07 (single 5-fold) and ~0.12 (single 70/30 holdout). Leaky = 1.00, fails by 0.30.
- B: 0.76 -0.15/+0.13 -> [0.61, 0.89]. Single 5-fold grouped: mean 0.72, sd 0.06, min 0.60 over 40 seeds. Leaky = 1.00 (fails by 0.11); chance 0.5 (fails by 0.11).
- Upper edge is below every leaky variant. Lower edge is above chance.

## Solver results (frozen)
correct (10x grouped 5-fold, selection inside): A 0.47, B 0.73 (truth 0.50 / 0.76; errors -0.03, -0.03; both pass). Variants that pass: 10-fold, 3-fold, leave-one-out/patient (A 0.44, B 0.76), single 5-fold C=0.1, C=10 unscaled, "A=0.5 by reasoning, B honest". Traps that fail: select-first (A only fails), random split (B only fails), both, chance 0.5/0.5 (fails B), selection-outside-with-grouped-folds (fails B), resubstitution.

## Known weaknesses / gaming / ambiguity
- **Chance guess** passes A by design; only B discriminates it. A is therefore weak evidence; it exists to catch selection-before-CV. Requires both fields to pass.
- **Single hold-out** (one 70/30 grouped split): ~32% fail B by noise (SD 0.125). Legitimate but noisy; not a bug in the grader.
- A model-free analyst could reason "A is noise so 0.5" and pass A without computing.
- Honest CV slightly underestimates (-0.03) because folds train on fewer patients; inside the band.
- Truth for B is conditional on this training set (the pipeline "as fitted on these data"); wording says "expected ROC AUC ... on patients not used for any fitting or selection", which an analyst reasonably reads as CV. Unconditional expectation is the same to 0.01.
- Feature-scale friction removed by stating values are normalised log-scale; features all have SD ~1. "z-score" and "C=1" are stated, but results are insensitive (checked C 0.1/10, no z-score).
- Rows in B are shuffled so unshuffled fold splitting does not accidentally group patients. If an agent uses `GroupKFold` it is correct; `KFold(shuffle=False)` on shuffled rows is a leak, as intended.
- `relevant_skills` picked by directory name only (no skill text read): `bio-machine-learning-model-validation`.
- Results depend on scikit-learn 1.9 / numpy 2.5; truth.json is recomputed by generate.py if versions change (float drift could shift the 4th decimal only; graded value is rounded to 2 dp).
