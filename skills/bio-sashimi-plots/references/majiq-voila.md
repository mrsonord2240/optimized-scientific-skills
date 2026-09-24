## MAJIQ-VOILA Interactive Viewer

**Goal:** Browse LSV posterior PSI distributions interactively with splice-graph topology.

**Approach:** Run `voila view` on MAJIQ output; it starts a local web server (open the printed address in a browser; there is no `-o` output file).

MAJIQ/VOILA (bundled with MAJIQ, majiq.biociphers.org) is licence-gated (academic/commercial download) and was **not installed or run** in testing; the commands follow MAJIQ's public docs, so check `voila view --help` for your version.

```bash
# MAJIQ V2 build: splicegraph.sql + one .voila file per quantification
voila view -p 5000 -j 8 build/splicegraph.sql psi_output/sample.psi.voila
voila view -p 5000 -j 8 build/splicegraph.sql deltapsi_output/group1_group2.deltapsi.voila

# MAJIQ V3 build (per the V2-to-V3 migration page): sg.zarr + the .psicov quantification + the .sgc coverage file of the group
voila view build/sg.zarr psi_output/Brain_Cerebellum.psicov build/Brain_Cerebellum.sgc
```

Do not mix V2 and V3 inputs in one call. Check `voila view --help` for the options (`-p`, `-j`) your version accepts.

VOILA shows:
- Complete LSV graphs (single source / single target nodes)
- Per-junction posterior PSI violin plots
- ΔPSI distributions across all conditions
- Confidence by junction within an LSV

**The only tool that visualizes complex multi-junction LSVs intuitively.** For events that don't fit canonical SE/A5SS/A3SS, VOILA is the visualization of choice. It needs the MAJIQ build's splicegraph plus the quantification file (`.voila` in V2; `.psicov` and `.sgc` in V3); without a licence use ggsashimi on the region instead.
