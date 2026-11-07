# Group-sequential and blinded SSR examples
```bash
Rscript "$SKILL_DIR/scripts/sequential.R" > sequential-example.txt
```
The entry point is scripts/sequential.R. Before running, set SKILL_DIR to the absolute skill directory and change into a newly created, empty working/output directory. The shell redirect truncates an existing sequential-example.txt, and the R plotting device can overwrite Rplots.pdf in the working directory.

**Done when:** the example output records group-sequential boundaries, survival design and blinded nuisance-only sample-size recalculation under the script's stated assumptions. This is an illustrative design, not trial-specific validation; review and adapt all parameters before use. Source versions target gsDesign 3.6+ and rpact 4.2+; check installed APIs. The survival example assumes proportional hazards; delayed effects require a suitable non-PH method or simulation. Blinded recalculation preserves error only under its stated conditions. See [methods](../references/group-sequential-ssr-enrichment.md).
