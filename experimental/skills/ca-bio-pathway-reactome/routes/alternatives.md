# Ranked and comparative workflows

```bash
Rscript examples/reactome_gsea.R
```

This source example uses a synthetic ranking and local human annotations, then writes temporary results and plots. It is not a user-data ranked CLI. For actual data adapt the complete signed, unique, decreasing Entrez ranking and fixed seed; inspect gene mapping, leading edges and correlation assumptions. `viewPathway` uses a pathway name for a local reaction view, while the web browser uses the stable ID. ReactomeGSA comparative/multi-omics methods contact a hosted service and require input-specific adaptation and authorization to transmit data. The source ORA demonstration is also synthetic. No individual diagnostic or treatment interpretation is supported.

Done when the method, ranking or design, database release, hierarchy overlap and any remote data transfer are documented.
