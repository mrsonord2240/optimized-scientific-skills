# scVI and scANVI
~~~sh
python scripts/scvi_scanvi_integration.py merged.h5ad integrated.h5ad --batch-key batch --max-epochs 100 --seed 0
~~~
Input X must be raw counts. Optional --labels-key cell_type enables scANVI and requires an Unknown label. Source defaults are 2000 HVGs, 10 latent dimensions, ZINB and CPU; output embeddings are coordinates, not corrected measurements or independently interpretable axes. Retain original counts outside the HVG-subset output. Check version compatibility and [model notes](../references/scvi-scanvi.md).
**Done when:** latent output, training settings, mixing and biological preservation are reviewed.
