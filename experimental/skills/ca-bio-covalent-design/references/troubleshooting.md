## Per-Tool Failure Modes

### Wrong warhead for residue

**Trigger:** A warhead/residue pairing is assumed from a broad class label.

**Mechanism:** Reaction depends on the residue microenvironment, electrophile, binding pose, and catalytic assistance; a class label does not establish residue selectivity.

**Symptom:** No covalent adduct observed despite docking pose.

**Fix:** Use literature-supported residue/warhead chemistry as a hypothesis, then verify site-specific adduct formation and competing reactivity experimentally.

### Excessive reactivity (off-target)

**Trigger:** Chloroacetamide in drug-candidate context.

**Mechanism:** Excess intrinsic electrophile reactivity can increase reaction with GSH and off-target nucleophiles.

**Symptom:** Toxicity in cell-based assays; non-specific binding signal.

**Fix:** Test a less intrinsically reactive electrophile and measure its GSH and target-reaction kinetics; alpha substitution can tune behavior but does not guarantee selectivity or stability.

### Geometric mismatch

**Trigger:** The ligand's reactive atom is poorly positioned relative to the target nucleophile.

**Mechanism:** Covalent reaction requires warhead-specific distance and approach geometry between the electrophilic atom and the nucleophilic atom (Cys Sγ for cysteine).

**Symptom:** No covalent labeling in mass spec despite predicted docking.

**Fix:** Identify the reaction atoms, inspect the pre-reaction Sγ-to-electrophile distance and reaction-specific angles, and use a docking protocol parameterized for that reaction. Do not substitute Cβ distance for the reacting sulfur.

### Reversibility unintended

**Trigger:** Designed irreversible TCI but warhead is reversible.

**Mechanism:** Reversibility depends on the complete electrophile, adduct chemistry, protein environment, and assay timescale; class-level labels are only hypotheses.

**Symptom:** Activity wanes after substrate washout in cellular assays.

**Fix:** Select chemistry with demonstrated behavior in the intended context and verify reversibility by dilution, washout, intact-protein MS, or another suitable experiment.

### kinact/Ki conflation

**Trigger:** Optimizing for IC50 instead of kinact/Ki.

**Mechanism:** Compounds with same IC50 differ in covalent efficiency.

**Symptom:** Compounds with similar endpoint IC50 values show different time-dependent target engagement or pharmacodynamic duration.

**Fix:** Fit a mechanistically appropriate kinetic model. Use `kinact/Ki` for qualifying irreversible two-step systems; use equilibrium, residence-time, or washout measurements where appropriate for reversible covalent systems.

### DOCKovalent over-prediction

**Trigger:** Default DOCKovalent run.

**Mechanism:** Covalent constraint forces docking; many ligands "succeed" but are unrealistic.

**Symptom:** Many compounds pass docking; few label in vitro.

**Fix:** Review measured/validated reactivity, reaction-atom geometry (Cys Sγ for cysteine), non-covalent recognition, strain, and site-specific experimental labeling.

## Reconciliation: Irreversible vs Reversible Covalent

| Aspect | Irreversible | Reversible covalent |
|--------|--------------|---------------------|
| Examples | KRAS G12C (acrylamide), BTK (ibrutinib) | Boronate (bortezomib), aldehyde (calpain inhibitors) |
| Toxicity profile | Off-target Cys labeling potential | Off-target equilibrium |
| Resistance mechanism | Mutation of reactive Cys | Mutation reduces affinity |
| Design decision | Consider duration of target engagement, safety, exposure, and resistance | Consider equilibrium, residence time, and recovery after washout |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Warhead not matching SMARTS | Different stereochemistry or charged | Use canonicalized + neutral SMARTS |
| DOCKovalent rejects ligand | No suitable Cys in pocket | Re-check residue accessibility |
| GSH adduct dominates | Warhead too reactive | Use less reactive warhead; or alpha-substitute |
| Off-target labeling in cells | Promiscuous warhead | Iterate warhead reactivity vs selectivity |
| Docking pose but no labeling | Geometric mismatch | Distance check; rotamer search |
| Intended irreversible inhibitor shows recovery after washout | Adduct chemistry is reversible or covalent reaction is incomplete | Re-evaluate the mechanism and fit the appropriate kinetic model |
| HCovDock fails on PROTAC | Tool optimized for monomer covalent | Use specialized tools for bivalent |
