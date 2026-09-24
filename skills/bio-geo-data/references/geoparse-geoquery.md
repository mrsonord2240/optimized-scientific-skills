# GEOparse and GEOquery

## GEOparse vs GEOquery

| Aspect | GEOparse (Python) | GEOquery (R/Bioconductor) |
|---|---|---|
| Maturity | OK; some known supplementary-file fetch issues since ~2022 | Mature; Bioconductor-supported |
| Output | `GEOparse.GSE` object with `gsms`, `gpls`, `metadata` dicts | `ExpressionSet` or list per platform |
| Supplementary files | `gse.download_supplementary_files()` (sometimes flakey) | `getGEOSuppFiles(gse)` (more reliable) |
| Integration | Pandas DataFrames | Bioconductor ecosystem |
| When | Python-first pipelines | R-first / use ExpressionSet downstream |

For production GEO workflows in R, GEOquery is the stable choice. For Python, GEOparse is the only option but verify file counts after download.

### GEOparse: full Series download

```python
import GEOparse


def get_gse(gse_id, dest='./geo_cache'):
    gse = GEOparse.get_GEO(geo=gse_id, destdir=dest)
    print(f'{gse_id}: {len(gse.gsms)} samples, {len(gse.gpls)} platforms')
    for gsm_name, gsm in list(gse.gsms.items())[:3]:
        print(f'  {gsm_name}: {gsm.metadata.get("title", ["?"])[0]}')
    return gse


# Supplementary files (raw data) -- verify file count manually after
gse = get_gse('GSE123456')
gse.download_supplementary_files(directory='./geo_cache')
```

### R: GEOquery (more reliable supplementary download)

```r
# Reference: Bioconductor GEOquery 2.70+ | Verify API if version differs
library(GEOquery)

gse <- getGEO('GSE123456', GSEMatrix = TRUE)
length(gse)             # one ExpressionSet per platform
head(pData(gse[[1]]))   # sample metadata
head(exprs(gse[[1]]))   # expression matrix (submitter-normalized -- verify processing notes)

# Raw / supplementary files
supp_dir <- getGEOSuppFiles('GSE123456', baseDir = './geo_cache')
list.files(rownames(supp_dir))
```
