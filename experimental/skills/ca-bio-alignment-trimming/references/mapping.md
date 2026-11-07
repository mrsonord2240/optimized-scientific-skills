# Mapping



When trimming for phylogenetic input, retain the column-index mapping so downstream site-specific analyses (selection per site, structure-mapped residues, etc.) can be back-traced:

```bash
trimal -in input.fasta -out trimmed.fasta -gappyout -colnumbering > kept_columns.txt

clipkit input.fasta -m smart-gap --log -o trimmed.fasta
```

These two flags emit different formats (checked on trimAl 1.4.1 and ClipKIT 2.14.0). `trimal -colnumbering` prints one line on stdout: `#ColumnsMap`, a tab, then a comma-separated list of the 0-based original column indices that survived (`#ColumnsMap	0, 1, 2, 4, ...`). `clipkit --log` writes `<output>.log` with one space-separated row per original column: 1-based position, `keep`/`trim`, site class (`constant`, `parsimony-informative`, `singleton`, `other`), gap fraction (e.g. `1 keep constant 0.0`). Parse both into the same 0-based kept-column list:

```python
kept_trimal = [int(x) for x in open('kept_columns.txt').read().split('\t', 1)[1].split(',')]
kept_clipkit = [int(row.split()[0]) - 1 for row in open('trimmed.fasta.log') if row.split()[1] == 'keep']
```
