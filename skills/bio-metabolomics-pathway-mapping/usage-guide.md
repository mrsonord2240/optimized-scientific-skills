# Pathway Mapping Usage Guide

## Overview

Pathway mapping places metabolomics results in biochemical context through over-representation (ORA), metabolite-set enrichment (MSEA/QEA), mummichog/PSEA on raw m/z features, and network-diffusion enrichment (FELLA). The central failure it guards against: enrichment destroys annotation uncertainty (a tentative ID becomes a confident p-value), the background set silently controls every result, topology "impact" is a hub artifact, and a steady-state pool size is not flux. The honest output is a hypothesis about network activity, conditional on the chosen annotations, background, database boundary, and ionization settings.

## Prerequisites

Install steps, the required local-session setup and the remote-call disclosures are in SKILL.md's Version Compatibility section. Know before starting: whether your metabolites are identified (IDs -> ORA/MSEA) or raw m/z features (-> mummichog/PSEA), the ionization mode and ppm of the run, and whether the compound list may leave the machine.

## Quick Start

Tell your AI agent what you want to do:
- "Run ORA on my list of identified significant metabolites against KEGG human pathways with an assay-specific background"
- "Predict pathway activity from my untargeted LC-MS feature table using mummichog with the full feature table as background"
- "Find the enzymes and reactions linking my metabolites using FELLA network diffusion"
- "Tell me whether my pathway result is being driven by a single hub metabolite"

## Example Prompts

### Identified-Compound Enrichment
> "Test which KEGG pathways are over-represented in my 18 confidently identified metabolites, using only the ~300 compounds my assay can detect as the background."
> "Run quantitative MSEA on my ranked metabolite fold changes instead of a cutoff-based ORA."

### Raw-Feature Activity Prediction
> "I have an untargeted negative-mode LC-MS peak table with m/z, p-value, and t-score but no IDs; predict perturbed pathways with mummichog and make sure the full table is the background."
> "Run integrated mummichog + GSEA PSEA at 5 ppm and report the predicted-active pathways as activity, not metabolite identifications."

### Mechanism and Sanity-Checking
> "Use FELLA diffusion to return the intermediate enzymes and reactions linking my KEGG compounds, and list which compounds did not map."
> "Check whether my top pathway is significant only because of L-alanine centrality, and re-run the enrichment across KEGG and SMPDB to see if it is robust."

## Related Skills

- metabolomics/metabolite-annotation - Annotation confidence levels (MSI) that feed ORA/MSEA and set the interpretive ceiling
- metabolomics/statistical-analysis - Upstream differential testing that produces the significant compound or feature list
- pathway-analysis/go-enrichment - Gene-set over-representation concepts
- pathway-analysis/gsea - Ranked-list enrichment concepts
- multi-omics-integration/mofa-integration - Joint gene+metabolite integration and its coverage-asymmetry traps
