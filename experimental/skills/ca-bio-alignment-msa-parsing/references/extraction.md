# Extract rows and regions

```python
from Bio import AlignIO
import sys
sys.path.insert(0, 'examples')
from msa_utils import normalize_alignment, select_columns

alignment = AlignIO.read('alignment.fasta', 'fasta')
ids = [r.id for r in alignment]
sequences = [str(r.seq) for r in alignment]
target = next((r for r in alignment if r.id == 'requested_id'), None)
region = alignment[0:10, 100:200]
ref = str(normalize_alignment(alignment)[0].seq)
ungapped_reference_columns = select_columns(alignment, [i for i, c in enumerate(ref) if c != '-'])
```

Choose the input, ID, row and column bounds and reference row to match the task. Coordinates are zero-based and end-exclusive. `select_columns` preserves aligned record and column annotations.
