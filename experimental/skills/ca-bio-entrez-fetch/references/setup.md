# Entrez setup for reference snippets

```python
from Bio import Entrez, SeqIO
Entrez.email = 'researcher@example.org'

def expect_start(text, prefix):
    if not text.lstrip().startswith(prefix):
        raise RuntimeError(f'Unexpected response: {text[:80]!r}')
    return text
```

Replace the email with the requester's address. Set `Entrez.api_key` only when a real key is available. Check response structure before parsing.
