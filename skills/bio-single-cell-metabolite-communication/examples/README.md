# Smoke-test fixture

Run python make_synthetic_fixture.py to create a deterministic 24-cell,
6-gene synthetic_metabolite_communication.h5ad. It contains the required
cell_type labels and uppercase human gene symbols, with a deliberately simple
Tumor enzyme / TCell sensor contrast.

Copy mebocost.conf.example to mebocost.conf and replace every placeholder
with the paths from the installed MEBOCOST database download. Then run the
guarded run_mebocost() pattern from SKILL.md against the generated h5ad.

Expected smoke-test invariants:

- AnnData has 24 observations, 6 variables, and a cell_type column.
- Tumor expresses higher NT5E and PTGES; TCell expresses higher PTGER2 and PTGER4.
- The script must be invoked through a __main__ guard on Windows.

Database versions and MEBOCOST scoring rules can change which communications
are significant, so this fixture deliberately does not promise a fixed FDR or
biological result. Confirm any result with the validation workflow in SKILL.md.
