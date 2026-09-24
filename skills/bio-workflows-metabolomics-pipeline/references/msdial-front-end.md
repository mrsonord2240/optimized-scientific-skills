# Alternative Front End -- MS-DIAL

When peak detection happens in the MS-DIAL GUI/console (MS2Dec deconvolution, GC-EI, DIA/SWATH), import the alignment-result table and enter the pipeline at Stage 2. The framing is unchanged: the imported table is still a parameterized hypothesis. See metabolomics/msdial-preprocessing for the export-parsing details, then continue with normalization-qc onward.

`scripts/msdial_import.R` turns the export into the objects Stage 2 expects (`feat`, `defs`, `sample_class`, `injection_order`, `batch_id`); run Stage 2 from `fm <- feat` on. Checked on a real MSDIALCUI 5.5.260820 export, whose 4 header rows (`Class`, `File type`, `Injection order`, `Batch ID`) precede the column header and whose `Class` cell marks where the per-sample columns start.

```r
source('scripts/msdial_import.R')   # export_file -> feat, defs, sample_class, injection_order, batch_id
```
