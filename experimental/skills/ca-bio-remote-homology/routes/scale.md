# Sequence collection search

```bash
mmseqs easy-search query.fa target.fa results.m8 tmp -s 7.5 --num-iterations 3
```

Replace inputs and use a fresh temporary directory. The source describes `-s 7.5` as a sensitive setting and reports a 5.7 default on its checked version. For clustering it gives `mmseqs easy-cluster proteins.fa cluster tmp --min-seq-id 0.5 -c 0.8`. For DIAMOND, build a database with `diamond makedb --in target.fa -d target` then search using `diamond blastp -d target -q query.fa -o results.tsv --more-sensitive -e 1e-10 -p 16 --outfmt 6`; the source also describes `--ultra-sensitive` and `blastx --frameshift 15` for applicable inputs. Verify CLI support and record thresholds.

Done when results and method-specific scores are retained with input and database provenance.
