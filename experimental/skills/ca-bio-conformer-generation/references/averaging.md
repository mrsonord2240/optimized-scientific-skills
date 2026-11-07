## Boltzmann Averaging of Properties

For ensemble descriptors (3D shape, dipole moment, polar surface area in 3D), Boltzmann-weight by energy:

```python
import numpy as np

def boltzmann_weights(energies, T=300.0):
    energies = np.array(energies)
    if energies.size == 0:
        raise ValueError('Boltzmann weighting requires at least one energy')
    if not np.all(np.isfinite(energies)):
        raise ValueError('Boltzmann weighting requires finite energies')
    if T <= 0:
        raise ValueError('Temperature must be positive')
    kt = 0.001987 * T  # kcal/mol at 300K
    rel = energies - energies.min()
    w = np.exp(-rel / kt)
    return w / w.sum()

def boltzmann_average(values, energies, T=300.0):
    if len(values) != len(energies):
        raise ValueError('values and energies must have the same length')
    w = boltzmann_weights(energies, T)
    return float(np.sum(np.array(values) * w))
```

Boltzmann populations require energies that approximate the relevant thermodynamic state, including conformer degeneracy and environmental effects where important. Raw gas-phase MMFF/UFF minima are exploratory surrogates, not generally validated populations.
