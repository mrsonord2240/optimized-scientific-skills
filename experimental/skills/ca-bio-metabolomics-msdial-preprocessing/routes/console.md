# Console preprocessing

```bash
MSDIALCUI.exe lcms -i input_folder -o output_folder -m method.txt
```

Replace paths and verify installed console help and data format. The source describes one `lcms` command for DDA and DIA, with acquisition type specified per file in an input CSV when needed; its DIA/ABF path was not validated on real DIA data. The `gcms` subcommand is the source's GC-EI route. Method files use `Key: Value`; a misspelled key may silently fall back to defaults. Archive parameters, acquisition metadata and outputs and verify the installed version's behavior before interpreting results. Console commands process data; the synthetic `examples/process_msdial_output.R` only demonstrates parsing and cleans its temporary files, not user data. For targeted concentration workflows use `bio-metabolomics-targeted-analysis`.

Done when acquisition mode, method, feature table, deconvolution and processing provenance have been independently checked.
