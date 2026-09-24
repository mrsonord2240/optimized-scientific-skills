# Koina-served predicted libraries

Moved from SKILL.md. Read when generating fragment intensities and iRT for a peptide list with Prosit/AlphaPeptDeep/MS2PIP/UniSpec through Koina.

## Generate a Predicted Library via Koina

**Goal:** Produce fragment intensities and iRT for a peptide list without a local GPU or wet-lab library.

**Approach:** Send a DataFrame of peptide sequences, charges, and collision energies to a Koina-hosted model; the dead proteomicsdb endpoint is replaced by the Koina server. Network calls are shown; the runnable example operates on an in-memory table so it needs no network.

```python
# Koina serves Prosit/AlphaPeptDeep/MS2PIP/UniSpec predictions; verify the
# koinapy constructor signature and input column names at runtime with help(Koina).
from koinapy import Koina
import pandas as pd

inputs = pd.DataFrame({
    'peptide_sequences': ['LGGNEQVTR', 'VEATFGVDESNAK'],
    'precursor_charges': [2, 2],
    'collision_energies': [30, 30]  # NCE; scan candidates and pick max spectral contrast
})

intensity_model = Koina('Prosit_2019_intensity', 'koina.wilhelmlab.org:443')
fragments = intensity_model.predict(inputs)  # mz, intensities, annotation per fragment

irt_model = Koina('Prosit_2019_irt', 'koina.wilhelmlab.org:443')
irt = irt_model.predict(inputs[['peptide_sequences']])  # arbitrary iRT units -- calibrate before use
```

Validate peptides before submitting: Prosit accepts only the 20 standard residues, length 7-30
(modified residues use the model's own notation, e.g. `M[UNIMOD:35]`). Anything else raises an
uncaught `tritonclient.utils.InferenceServerException` from the server rather than a clear message.

```python
import re
from tritonclient.utils import InferenceServerException

def valid_prosit_peptide(seq, lo=7, hi=30):
    return lo <= len(seq) <= hi and re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWY]+', seq) is not None

bad = inputs[~inputs['peptide_sequences'].map(valid_prosit_peptide)]  # report and drop these
try:
    fragments = intensity_model.predict(inputs.drop(bad.index))
except InferenceServerException as e:
    raise RuntimeError(f'Koina rejected the request: {e}') from e
```
