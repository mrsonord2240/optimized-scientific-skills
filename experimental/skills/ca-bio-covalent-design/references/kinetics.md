## Kinetics: kinact / Ki

Covalent inhibition kinetics:
- **Ki**: reversible binding affinity (initial, like non-covalent IC50)
- **kinact**: rate of covalent bond formation (sec^-1)
- **kinact/Ki**: second-order rate constant, "covalent efficiency" (M^-1 s^-1)

For irreversible two-step inhibition that follows the corresponding kinetic model, report fitted `kinact`, `Ki`, and `kinact/Ki` rather than only a time-dependent IC50. Two compounds with the same IC50 can have different kinetic components:
- Low Ki, low kinact: tight binding, slow covalent bond
- High Ki, high kinact: loose binding, fast covalent bond

Because `kinact/Ki` depends on the target, construct, assay conditions, and kinetic model, compare values within a matched assay series and alongside exposure, intrinsic reactivity, target engagement, and selectivity. Reversible-covalent systems may require different mechanistic models and residence-time or washout measurements.

## Intrinsic Reactivity Assays

Before committing to a warhead, measure intrinsic reactivity (off-target risk):

```python
# Generic GSH stability assay readout - measure half-life of warhead with 10 mM GSH
# kinact_GSH from time-course of warhead disappearance
```

Do not assign a universal GSH half-life from the warhead name alone. Substitution, electronics, ionization, solubility, and assay conditions can change the observed rate. Report the GSH concentration, buffer, temperature, analytical method, and fitted half-life or second-order rate constant for the complete compound.
