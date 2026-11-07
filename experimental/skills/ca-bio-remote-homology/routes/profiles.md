# Iterative profiles

```bash
psiblast -query query.fa -db database -num_iterations 3 -inclusion_ethresh 0.002 -evalue 0.01 -num_threads 8 -out_pssm profile.asn -out_ascii_pssm profile.txt -out results.txt
```

Supply a formatted protein database and query. Save the PSSM for reuse with `-in_pssm`; inspect included hits for profile drift. Alternatively `jackhmmer -N 3 --chkhmm iter --tblout hits.tbl query.fa database.fa` saves checkpoints; early convergence can stop before the third round. The source notes 0.005 as PSI-BLAST's inclusion default and 0.002 as a stricter example. `examples/iterative_profile.sh query.fa database.fa` compares three methods; it builds local databases and writes fixed names, so use a fresh working directory and verify database inputs. Counts from different methods are not directly comparable.

Done when input, database version, iterations, scores and limitations are recorded and included hits reviewed.
