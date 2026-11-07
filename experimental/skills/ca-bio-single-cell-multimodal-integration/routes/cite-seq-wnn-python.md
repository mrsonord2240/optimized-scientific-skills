# CITE-seq in Python: CLR and Muon WNN

```bash
python scripts/cite_seq_muon.py filtered_feature_bc_matrix.h5 cite_out --resolution 0.5
```

A CLR-only fallback: QC, per-modality PCA and neighbours, Muon WNN, Leiden on the WNN graph, UMAP. It writes `cite_out.h5mu` and a UMAP PDF.

- CLR does not remove ADT background and Python has no first-party DSB. For denoised protein use `routes/cite-seq-wnn-r.md` (DSB) or `routes/totalvi.md`, or run DSB in R and re-import.
- Muon WNN needs modality-local graphs built first; the script does that. WNN lives on the `MuData` object: use `mu.tl.umap(mdata, neighbors_key='wnn')`, not `scanpy.tl.umap`.
- After per-modality QC the modalities hold different cells: `mu.pp.intersect_obs(mdata)` before any paired analysis (the script does). Editing `mdata.mod['rna'].obs` needs `mdata.update()` to reach `mdata.obs`.
- Muon's default of 200 candidate neighbours exceeds small datasets; the script caps it.

Done when `cite_out.h5mu` exists and the answer says protein was not denoised.
