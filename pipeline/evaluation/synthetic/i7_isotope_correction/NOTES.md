# i7_isotope_correction (not agent-visible)

## Pitfall and source
Raw mass-isotopologue distributions (MIDs) contain natural-abundance isotopes of every element in the measured ion
(13C incl. derivatisation carbons, 2H, 15N, 17O/18O, 29Si/30Si) and tracer impurity. Enrichment claims need full correction
(Buescher et al. 2015 Curr Opin Biotechnol; Midani et al. 2017). Design is from that pitfall, not from any Skill.

## World (generate.py, seed 20260919)
- Tracer [U-13C6]-glucose, 98 atom% 13C at each labelled position (98 chosen, not 99, so ignoring purity fails by margin at high enrichment).
- 4 metabolites; measured ion formula given: lactate 2TBDMS [M-57] C11H25O3Si2 (n=3), alanine 2TBDMS C11H26NO2Si2 (n=3),
  glutamate 3TBDMS C19H42NO4Si3 (n=5), citrate [M-H]- C6H7O7 (n=6). Si and 8-14 unlabellable derivative carbons make non-13C
  correction large; citrate (underivatised) is the low-offset case.
- 4 time points (1, 4, 8, 24 h) + one unlabelled 24 h control per metabolite. True corrected MIDs (number of tracer-derived carbons)
  are written explicitly in METS. Enrichments 0.02-0.76.
- Forward model: column k = Binom(k, purity) 13C at labelled positions (x) natural 13C on the other Ctot-k carbons (x) natural H, N, O, Si.
  Truncated at M+(n+2), random per-sample scale (2e5-5e7), 0.5% relative Gaussian noise per peak. Unmeasured columns = NA.
- Enrichment definition (stated in the task): sum_i i*MID_i / backbone_carbons, where i counts tracer-derived carbons (a tracer-derived
  carbon that is 12C because of impurity still counts). Same convention as isocor.
- Truth = constructed MIDs, rounded to 4 dp, written into eval.json by generate.py. Answer keys: `<metabolite>_<time>h`.

## Grader and tolerance
One numeric_tolerance, +/-0.01 absolute on each of 6 values. Reference correction error vs truth: max 0.0022 on the real data;
Monte-Carlo over 300 noise redraws x 6 values: mean 0.0007, p99 0.0027, max 0.0040. Variants also pass: older abundance table
(0.989/0.011 13C, Si 4.67/3.1%): max 0.003; square M0..Mn solve: max 0.002; isocor 2.2.4 (moiety + derivative_formula, purity [0.02,0.98]): max 0.0022.
So tolerance is ~2.5x the worst legitimate-variant error and about 1x the effect of the smallest trap.

## Results (`python solvers.py`), fields failing out of 6 (|err| > 0.01)
| solver | fail | note |
| --- | --- | --- |
| correct | 0 | |
| raw MID | 6 | errors 0.01-0.11; citrate only 0.012 (thin) |
| 13C-only correction (all ion carbons) | 4 | passes lactate_24h and citrate_8h |
| correction with underivatised metabolite formula | 5 | correct for citrate by construction |
| ignore tracer purity | 1 | lactate_24h err -0.018 (only failing field; 1.8x tol) |
| subtract control raw enrichment | 1 | lactate_24h err -0.034 (only failing field) |
| subtract control MID fractions | 4 | |
| n = carbons in ion formula | 5 | |

## Known weaknesses / blind spots
- The purity trap and the control-subtraction trap each fail on ONE field only (lactate_24h). They are still caught by the all-or-nothing
  grade, but a solver that gets lactate_24h right by luck would pass them. Purity is small at 98% unless enrichment is high.
- Underivatised citrate barely discriminates: raw error 0.012, 13C-only error 0.005 (passes). Kept so the requested set mixes regimes.
- The 13C-only trap only fails because of the Si/derivative atoms in three of the four ions. On a purely underivatised LC-MS set (H/N/O only)
  it would pass at +/-0.01 (error ~0.005): that is a real blind spot of the design if Si is absent.
- The task states the enrichment convention ("attributable to the tracer", "tracer-derived carbon", backbone_carbons) and gives
  tracer_13C_purity and the ion formula; it does not say to correct. A solver that treats measured 13C content as enrichment (no purity
  step) gets the "ignore purity" answer.
- Ambiguity risk: derivatisation moiety vs metabolite moiety split (isocor wants them separately). The ion formula is given as a single
  formula; the analyst must work out that backbone_carbons < C in the formula. Formula strings were verified to add up (moiety + derivative = ion).
- Gaming: raw M+n peak patterns could be fit without the formula by using the control sample as an empirical natural-abundance pattern;
  a correct solver using the control's M0..Mn pattern as the unlabelled column would still be approximately right for k=0 only, so not exact.
  Not tested.
- MIDs are hand-set (not from a kinetic model); realistic in shape, not in dynamics.
