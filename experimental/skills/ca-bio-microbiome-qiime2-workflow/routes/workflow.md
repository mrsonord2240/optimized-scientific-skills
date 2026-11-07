# Orchestrate a study

```bash
qiime metadata tabulate --m-input-file metadata.tsv --o-visualization metadata.qzv
```

Inspect IDs and metadata types before wiring actions. Use a `#q2:types` row when integer-like IDs need categorical treatment. The source `examples/qiime2_16s.sh` is a fixed-path orchestration skeleton: it imports, denoises, classifies, builds a tree, runs diversity and ANCOM-BC, replays and exports. Its truncation and sampling-depth values are placeholders; do not run unchanged on a study. Adapt classifier and scientific parameters after checking the installed release and inputs. Classifier compatibility depends on its training environment. Defer scientific decisions to `bio-microbiome-amplicon-processing`, `bio-microbiome-taxonomy-assignment`, `bio-microbiome-diversity-analysis`, `bio-microbiome-differential-abundance` and `bio-microbiome-functional-prediction`. Outputs can collide with existing files or directories.

Done when metadata, parameters, artifact types, exclusions, tool versions and provenance are reviewed for the study.
