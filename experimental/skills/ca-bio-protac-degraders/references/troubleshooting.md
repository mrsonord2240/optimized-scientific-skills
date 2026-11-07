## Per-Tool Failure Modes

### PRosettaC -- inaccessible E3 in selected ligase

**Trigger:** Target's known binding mode incompatible with E3 ligase orientation.

**Mechanism:** The selected binary poses, exit vectors, linker conformations, or protein orientation may not support a compatible ternary geometry.

**Symptom:** Low ternary complex scores; high RMSD across replicates.

**Fix:** Try a different E3 such as CRBN or VHL and compare against experimentally resolved ternary complexes with compatible exit-vector geometry.

### DeepTernary -- novel chemotype

**Trigger:** Target ligand or E3 ligand outside training distribution.

**Mechanism:** A ligand, linker, target, or E3 outside the method's validated domain may require extrapolation.

**Symptom:** Predicted ternary complex unrealistic.

**Fix:** Compare with an independently configured structural method and relevant known complexes; validate prospective ranking experimentally.

### Hook effect at high PROTAC concentration

**Trigger:** PROTAC concentration becomes high enough that separate target-PROTAC and E3-PROTAC binary complexes compete with productive ternary-complex formation.

**Mechanism:** Saturation by binary complexes reduces the population of productive ternary complex. Negative cooperativity can worsen ternary formation but is not the definition of the hook effect.

**Symptom:** Degradation increases and then decreases across a sufficiently broad concentration-response experiment.

**Fix:** Confirm the downturn experimentally over a broad dose range, then optimize ternary-complex geometry, cooperativity, exposure, and dosing without assuming that linker shortening alone will solve it.

### Insufficient cell permeability

**Trigger:** Measured permeability or cellular exposure is poor relative to biochemical activity.

**Mechanism:** Size, exposed polarity, conformation, ionization, or efflux may limit intracellular exposure.

**Symptom:** Cellular degradation potency is substantially worse than biochemical ternary-complex or binding measurements.

**Fix:** Optimize linker and exposed polarity using measured permeability, solubility, and intracellular exposure across the series. Do not impose a universal MW or TPSA cutoff.

### E3-target distance miscalculation

**Trigger:** Computing linker length from binary models without ternary refinement.

**Mechanism:** A distance from separately positioned binary structures does not determine the accessible ternary geometry or linker conformational ensemble.

**Symptom:** PROTACs synthesized at wrong linker length; no degradation.

**Fix:** Use a ternary structural hypothesis to define a chemically diverse linker series, then compare conformational feasibility and experimental degradation across that series.

### Molecular glue vs PROTAC confusion

**Trigger:** Designing as PROTAC when target lacks defined ligand.

**Mechanism:** A molecular glue stabilizes or induces a protein-protein interaction without the two-ligand-plus-linker architecture assumed by a PROTAC workflow.

**Symptom:** Design too rigid; no degradation despite ternary prediction.

**Fix:** For targets without known ligand, consider molecular glue discovery instead.


## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| PRosettaC fails to converge | Input geometry, sampling, or service/configuration problem | Inspect inputs and logs; compare justified recruiter/linker hypotheses |
| DeepTernary returns clashing pose | Prediction outside a validated domain or incorrect interface | Inspect confidence and clashes; compare an independent method and known structures |
| AlphaFold3 ternary unrealistic | Unrestrained prediction has a low-confidence or incorrect interface | Inspect confidence and compare with PRosettaC or known ternary structures |
| Cellular phenotype disagrees with target-degradation assays | Exposure, off-target degradation, assay timing, or pathway effects | Measure target engagement/degradation and use proteome-wide selectivity assays where appropriate |
| Degradation decreases at high PROTAC concentration | Hook effect from competing binary complexes | Confirm with a broad dose range; optimize ternary geometry and exposure |
| Synthesis is impractical | Proposed connectivity lacks a credible route | Obtain medicinal-chemistry review and redesign attachment chemistry or linker |
| Poor permeability or intracellular exposure | Size, exposed polarity, conformation, or efflux | Measure the bottleneck and optimize the series; avoid a universal size cutoff |
| `protac_enumerate.build_protac()` raises "Attachment dummy bonds must be single bonds" | An explicit `[*]` exit-vector dummy is attached via a double bond or other non-single bond | Use a single-bond attachment point; only relax this if the intended connection chemistry has its own, separately validated bond-order rule |
