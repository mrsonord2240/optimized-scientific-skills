# Randomize two treatments within blocks
~~~sh
Rscript scripts/randomize_blocks.R units.csv assignment.csv ctrl,treat 20260528
~~~
Input CSV needs unique id and block columns. The script randomizes the supplied two labels within each block, then independently randomizes order within each block. Even block sizes give equal counts; odd sizes have one extra of the first label, so review allocation and feasibility in advance. The final argument is the recorded seed. Constraints beyond this simple blocked scheme require a study-specific allocator.
**Done when:** unit counts, per-block balance, processing order, constraints and seed have been reviewed and saved.
