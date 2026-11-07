---
name: ca-bio-virtual-screening
description: Use when preparing or screening molecular ligands against a protein structure by docking.
category: Data Analysis
license: MIT
author: GPTomics
---

# Virtual Screening

Choose a route, read it, and run its command before writing new code. Paths are relative to this Skill.

| Need | Read |
|---|---|
| Prepare receptor or ligand | `routes/preparation.md` |
| Identify a binding pocket | `routes/site.md` |
| Dock a single ligand | `routes/docking.md` |
| Screen a ligand library | `routes/library.md` |
| GNINA scoring or large screening | `routes/large-screen.md` |
| Convert docked PDBQT poses to SDF | `routes/poses.md` |
| Related analysis | `routes/related.md` |

Choose and validate receptor, pocket, protonation, search box and controls; docking scores and poses require independent validation.
Source checked: Vina 1.2.5+, meeko 0.5+, RDKit 2024.09+.
