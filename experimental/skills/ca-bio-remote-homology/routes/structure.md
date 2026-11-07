# Structure comparison

```bash
foldseek easy-search query.pdb structure_db results.m8 tmp
```

Supply existing structure and compatible database inputs. The source also describes ProstT5 for sequence-only 3Di search and Foldseek multimer modes; check installed help and model requirements before adapting. `examples/foldseek_search.sh query.pdb db_dir tmp_dir` may download databases and writes `results.m8`; inspect storage and use a fresh working directory. Structural similarity alone does not establish homology. Report actual scores and uncertainty alongside supporting evidence.

Done when database release, alignment statistics, coverage and interpretation limits are recorded.
