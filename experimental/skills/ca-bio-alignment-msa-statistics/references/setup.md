# Set up source metric functions

```python
import sys
sys.path.insert(0, 'examples')
from msa_utils import load_alignment, check_alphabet, pick_background
alignment = load_alignment('alignment.fasta')
background, label = pick_background(alignment)
check_alphabet(alignment, background)
```

Select the actual input file and verify its inferred format and alphabet. The loader normalizes case, `.` and `~` gaps, and maps RNA U to T for nucleotide data. Preserve case separately for A2M/A3M insert states.
