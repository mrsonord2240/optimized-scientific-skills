# Environment compatibility
Verified 2026-09-30. No single environment satisfies every tool; keep three:

| Env | Packages | Used for |
|---|---|---|
| `chrombpnet` (Python 3.8) | chrombpnet 1.0.1 (pins tensorflow 2.8.0, keras 2.8.0, numpy 1.23.4), kundajelab/variant-scorer (git 0e1e341), pybedtools (conda; the pip build fails), bedtools and UCSC `bedGraphToBigWig` (conda `bedtools`, `ucsc-bedgraphtobigwig`), MEME suite (`tomtom`) | training, QC, `pred_bw`, `contribs_bw`, variant-scorer |
| `torch` (Python 3.11) | torch >= 2.1 (cu128 build for RTX 50-series), bpnet-lite 1.0.0, tangermeme 1.5.0, modisco-lite 2.4.0, enformer-pytorch 0.8.12, pyfaidx | attributions, TF-MoDISco, variant effect, Enformer |
| `scbasset` (Python 3.11) | tensorflow-cpu 2.15.1, keras 2.15.0, numpy 1.26.4, pandas < 3, anndata 0.11.4, psutil, scBasset (git aed3a6f, not on PyPI) | scBasset only |

TensorFlow 2.8 has no CUDA 12 build, so chromBPNet trains on CPU (the GPU helps only in the `torch` env). The PyPI name of TF-MoDISco-lite is `modisco-lite`. Deep-learning tooling evolves rapidly; verify with `pip show <package>` and `<tool> -h`.
