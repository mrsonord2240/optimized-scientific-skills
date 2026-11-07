# VST, rlog and normTransform (visualization only)

## VST vs rlog vs normTransform (Visualization Only)

**Goal:** Produce homoskedastic log2-scale counts for PCA, heatmaps, clustering, ML features.

**Approach:** Use `vst()` by default. Switch to `rlog()` only for n<30 with size factors varying >4x. Never use `normTransform()` for distance-based plots. Never use VST/rlog values as input to DE.

```r
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
```

The `vst()` function default is `blind=TRUE`, but the current DESeq2 vignette recommends `blind=FALSE` for any downstream visualization AFTER the model is fit (it uses the design when fitting dispersions; appropriate when the design is already settled). Reserve `blind=TRUE` for unsupervised QC where the design should not influence the transformation (e.g., "is this sample consistent with its group?"). The vignette has flip-flopped over the years on which to recommend by default -- pass `blind=` explicitly.

`vst()` fits the dispersion trend on `nsub = 1000` genes chosen deterministically to span the range of mean normalized counts (not the most variable genes). With fewer than 1000 genes it errors ("less than 'nsub' rows"): set `nsub` lower or call `varianceStabilizingTransformation()` directly.
