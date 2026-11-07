# Choosing a method

Anchor structure first:

| Anchor structure | What is shared | Example | Methods | Route |
|------------------|----------------|---------|---------|-------|
| Vertical / paired | Same cells | CITE-seq, 10x Multiome | WNN, totalVI, MultiVI (paired), MOFA+ | `routes/cite-seq-wnn-r.md`, `routes/totalvi.md`, `routes/multiome-wnn.md`, `routes/mofa.md` |
| Diagonal / unpaired | Nothing; a prior feature graph | Independent scRNA + scATAC | GLUE, Seurat v5 bridge | `routes/glue.md`, `routes/bridge.md` |
| Mosaic | Some modalities only | Batch A RNA + ATAC, batch B RNA | MultiVI, StabMap, Cobolt | `routes/multivi.md` |
| Horizontal | Same modality, different cells | Batch correction | not this Skill | `routes/other-skills.md` |

Two paired datasets sharing only one modality (a Multiome and a CITE-seq run sharing RNA) are a mosaic problem: anchor on the shared RNA and impute or bridge the rest (MultiVI, StabMap, bridge integration), do not force one WNN.

Paired methods:

| Method | Use when | Fails when |
|--------|----------|------------|
| WNN (Seurat) | Fast joint embedding and clustering of one well-normalized paired dataset | Protein background not removed first; a noisy or saturating modality dominates; unpaired or mosaic data |
| totalVI | Denoised protein, DE, batch integration, merging antibody panels | Tiny datasets (VAE overfits); very large data without GPU |
| MultiVI | Paired plus unpaired RNA/ATAC (mosaic); generative DE/DA | `batch_key` set to the modality indicator; imputed modality read as measured |
| MOFA+ | Interpreting shared versus modality-specific variation | Used for clustering or denoising; likelihood mismatched to data |
| mojitoo | Quick paired joint reduction from existing PCA/LSI | No way to down-weight a noisy modality; paired only |

Unpaired and mosaic methods:

| Method | Use when | Fails when |
|--------|----------|------------|
| GLUE | Diagonal scRNA + scATAC; regulatory inference as a by-product | Genome-build mismatch (empty guidance graph); adversarial over-mixing of distinct states |
| Seurat v5 bridge | Map a query onto a reference in another modality via a multiome bridge | Poor or batch-mismatched bridge; rare query-only populations mislabelled |
| StabMap | Mosaic with informative unshared features | Error compounds per hop |
| Cobolt, scMoMaT | Mosaic with a generative shared latent | DE on imputed values |

ADT normalization: CLR (Seurat `NormalizeData(method='CLR', margin=2)`) is quick and needs no empty droplets but does not remove background. DSB needs the raw matrix with empty droplets and removes ambient and per-cell technical noise. The CLR margin convention differs across versions; check `?NormalizeData`. When methods compete, check the installed tool docs for the current default.

Done when one route is chosen.
