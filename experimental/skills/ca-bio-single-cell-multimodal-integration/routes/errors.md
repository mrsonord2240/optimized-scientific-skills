# Multimodal integration problems

| Symptom | Cause | Fix |
|---------|-------|-----|
| WNN clustering driven by ADT | A few saturating high-variance proteins dominate the graph | Report weights; denoise ADT (DSB or totalVI); recheck clustering stability |
| Background smear in every ADT cluster | CLR or WNN without empty-droplet denoising | `routes/cite-seq-wnn-r.md` or `routes/totalvi.md` |
| DSB nonsense output, no error | Filtered matrix passed as the empty droplets | Supply the raw matrix empties; the script checks the median |
| Spurious batch structure after merging Multiome runs | Per-dataset peak sets | Re-quantify all cells against one peak set |
| GLUE blob, no alignment | Near-empty guidance graph from a build or coordinate mismatch | One genome build for genes and peaks |
| RNA and protein disagree for a marker | Often real post-transcriptional biology | Treat as informative; do not correct away |
| MultiVI batch effects persist | `batch_key` was the modality indicator | A separate covariate for real batches |
| DE on a modality looks too clean | Run on imputed or gene-activity values | Flag as model-dependent; validate on a measured modality |
| `scanpy.tl.umap` fails on Muon WNN | Muon WNN metadata is MuData-specific | `mu.tl.umap(mdata, neighbors_key='wnn')` |
| Global `obs` stale after editing a modality | MuData does not propagate | `mdata.update()` |
| Cell sets differ between modalities | Per-modality QC | `mu.pp.intersect_obs(mdata)` first |
| R and Python object round-trip loses information | MuDataSeurat and zellkonverter are lossy | Stay in one ecosystem |

Done when the symptom's fix has been applied and the script rerun.
