# Local DALI with DaliLite v5

Local DALI with DaliLite v5 (the DALI web server is the alternative). PDB ids for `import.pl` must be 4 characters; chains become `1mbnA`, `1a3nA`, ...:

```bash
import.pl --pdbfile 1mbn.pdb --pdbid 1mbn --dat DAT/
import.pl --pdbfile 1a3n.pdb --pdbid 1a3n --dat DAT/
dali.pl --cd1 1mbnA --cd2 1a3nA --dat1 DAT/ --dat2 DAT/ --title mb_hb --outfmt summary
# results: ./1mbnA.txt  ->  "1:  1a3n-A  20.3  1.6  141  141  26  MOLECULE: HEMOGLOBIN (ALPHA CHAIN)"
#                                          Z   rmsd lali nres %id
```

The result is written to `<cd1>.txt` in the working directory (`--title` only names the job) and the next job with the same query chain overwrites it, so run each job in its own directory. An empty table means no significant similarity. Checked: myoglobin vs hemoglobin alpha Z 20.3 (26% id); PKA (1ATP chain E) vs CDK2 (1HCK) Z 24.1; myoglobin vs PKA empty.
