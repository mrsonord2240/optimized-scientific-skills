## Isobaric (TMT/iTRAQ) Quantification

### Extract and impurity-correct reporter ions

**Goal:** Pull TMT reporter intensities from spectra and correct cross-channel isotope bleed.

**Approach:** Read spectra on disk, `quantify` the reporter region, then `purityCorrect` with a LOT-SPECIFIC impurity matrix from the reagent Certificate of Analysis. `readMSnSet` reads an already-quantified text matrix and does NOT extract reporters.

```r
library(MSnbase)

raw <- readMSData('experiment.mzML', mode = 'onDisk')
# method='max' for centroided spectra; reporters=TMT10 defines the 126-131 reporter m/z.
# MSnbase 2.32.0 has no TMT18 set.
quant <- quantify(raw, reporters = TMT10, method = 'max')
# TMTpro 16plex (126..134N): swap in TMT16 for TMT10 here and use the CoA route below
# quant <- quantify(raw, reporters = TMT16, method = 'max')

# edit = FALSE: the default edit = TRUE calls edit(M) and blocks under Rscript / on a cluster.
# x = is a MANUFACTURER TEMPLATE and MSnbase ships templates only for x = 4, 6, 8, 10; any other x
# (TMTpro 16) falls through to an unnamed diag(x) and stops with "length of 'dimnames' [1] not equal
# to array extent". REPLACE with lot-specific Certificate of Analysis values -- for TMTpro the only route:
# imp <- makeImpuritiesMatrix(filename = 'lot_coa.csv', edit = FALSE)
#   CSV layout (MSnbase extdata TMT6plexPurityCorrections.csv): a leading Tag column of channel names
#   (read as row.names; omit it and R stops with "duplicate 'row.names' are not allowed"), then n
#   neighbour-OFFSET columns in percent (-n/2..-1, +1..+n/2), so a 16plex CoA has Tag + 16 columns:
#     Tag,-3,-2,-1,+1,+2,+3
#     126,0,0,0,6.1,0,0
#     127,0,0,0.5,6.7,0,0
imp <- makeImpuritiesMatrix(x = 10, edit = FALSE)
quant <- purityCorrect(quant, imp)
```

### Bridge multiple TMT plexes with IRS

**Goal:** Make reporter intensities comparable across separate TMT runs.

**Approach:** Absolute reporter intensities for the same protein differ 2-5x between plexes because each plex samples a random point on the elution profile. Sample-loading normalization fixes within-run loading; the Internal Reference Scaling bridge (Plubell 2017) then pins each plex's pooled reference channel to a common per-protein value. Order: SL, then IRS. A protein whose reference is 0 or missing in any plex cannot be bridged: mask and report it rather than letting it become Inf/NaN. Check the bridge on the NON-reference channels (IRS forces the reference channels equal by construction).

Code: `sample_loading_normalize` and `irs_scale` in `examples/lfq_normalization.py` (input: protein x channel table of summed PSM reporter ions, one pooled reference channel per plex given to `irs_scale` as a column name per plex; the example plants a plex offset and prints it before and after IRS). Use `plex_offset` there to check the bridge on the non-reference channels.

Worked example with a planted plex offset: `examples/lfq_normalization.py`.
