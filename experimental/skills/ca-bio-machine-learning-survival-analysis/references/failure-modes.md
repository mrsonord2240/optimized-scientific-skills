# Survival Modeling Failure Modes

Survival-model failure modes: C-index-only reporting, Kaplan-Meier under competing risks, Fine-Gray misreading, immortal time, selection-before-CV.

## Per-Method Failure Modes

### Reporting only the C-index
- **Trigger:** Summarizing a survival model by Harrell's C alone.
- **Mechanism:** C is censoring-dependent, monotone-invariant (blind to calibration), and insensitive.
- **Symptom:** Great C, badly miscalibrated absolute risks; ranking-equivalent models look identical.
- **Fix:** Uno's C(tau) + AUC(t) + IBS-vs-KM + calibration curves on out-of-sample data.

### Kaplan-Meier under competing risks
- **Trigger:** Using 1-KM for incidence when a competing event exists.
- **Mechanism:** 1-KM assumes competing-event subjects could still have the event; they cannot.
- **Symptom:** Incidence overestimated; sums across causes exceed 1.
- **Fix:** Report the CIF (Aalen-Johansen); use cause-specific and Fine-Gray models.

### Misverbalizing a Fine-Gray coefficient
- **Trigger:** Saying a Fine-Gray sHR "increases the rate of the event."
- **Mechanism:** The subdistribution hazard maps to cumulative incidence, not the event rate.
- **Symptom:** Causal/etiologic claims from a prognostic model.
- **Fix:** State it as an effect on *cumulative incidence*; report cause-specific too.

### Immortal time bias
- **Trigger:** Grouping subjects by a post-baseline event.
- **Mechanism:** Guaranteed event-free time is attributed to the exposed group.
- **Symptom:** Spectacular development performance, external collapse.
- **Fix:** Time-varying exposure, landmarking, or target-trial emulation.

### Selection-before-CV in a signature
- **Trigger:** Selecting genes on all data, then CV-ing the final Cox model.
- **Mechanism:** Held-out folds informed selection (the dominant overfitting capacity in p>>n).
- **Symptom:** Optimistic, irreproducible signature.
- **Fix:** Selection inside nested CV with survival metrics; elastic net for stability.
