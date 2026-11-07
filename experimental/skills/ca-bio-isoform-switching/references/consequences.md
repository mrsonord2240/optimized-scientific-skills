## Functional Consequence Annotation

**Goal:** Predict how each switch alters protein structure, function, and stability.

**Approach:** Get ORFs, extract sequences, run the external annotators outside R, import their results, then request only the consequence types whose annotation was imported.

Order (each step needs the previous one; verified on a synthetic GTF without CDS and on the real chrX GTF with CDS):

1. **ORFs first.** `importRdata(addAnnotatedORFs = TRUE)` loads CDS from a GTF that has CDS lines (`orf_origin` = Annotation). Only if the GTF has no CDS (or isoforms are novel) run `analyzeORF(aSwitchList, orfMethod = 'longest', genomeObject = NULL)`. `extractSequence()` and `analyzePFAM()` stop without ORFs ("Please run the 'addORFfromGTF()' ... function(s) to detect ORFs").
2. **Do not run `analyzeORF()` over annotated ORFs** (for novel isoforms the package vignette adds `addORFfromGTF()` + `analyzeNovelIsoformORF()`; not run here). It overwrites all of them (real chrX: `orf_origin` Annotation 1738 -> Predicted 1738) and worsens the NMD call: PTC flag versus Ensembl `nonsense_mediated_decay` biotype was sensitivity 0.96 / specificity 1.00 with annotated ORFs and 0.64 / 0.95 after `analyzeORF('longest')`.
3. `analyzeAlternativeSplicing(aSwitchList, onlySwitchingGenes = TRUE)` — needed for `intron_retention`.
4. `extractSequence(aSwitchList, onlySwitchingGenes = TRUE, pathToOutput = 'sequences/', writeToFile = TRUE)` (create the directory first) writes `isoformSwitchAnalyzeR_isoform_nt.fasta` and `..._AA.fasta`.
5. Run the annotators (table) and import each result.
6. `analyzeSwitchConsequences()` with only the types you imported.

```r
# ORFs: annotated CDS come with importRdata; predict only when the GTF had none
if (!any(aSwitchList$orfAnalysis$orf_origin == 'Annotation', na.rm = TRUE)) {
    aSwitchList <- analyzeORF(aSwitchList, orfMethod = 'longest', genomeObject = NULL)
}
aSwitchList <- analyzeAlternativeSplicing(aSwitchList, onlySwitchingGenes = TRUE)
dir.create('sequences', showWarnings = FALSE)
aSwitchList <- extractSequence(aSwitchList, onlySwitchingGenes = TRUE, pathToOutput = 'sequences/', writeToFile = TRUE)

# after running the external tools on sequences/:
aSwitchList <- analyzeCPC2(aSwitchList, pathToCPC2resultFile = 'cpc2_result.txt', removeNoncodinORFs = FALSE)
aSwitchList <- analyzePFAM(aSwitchList, pathToPFAMresultFile = 'pfam_scanfmt.txt')
# licence-gated, not run here: analyzeSignalP(pathToSignalPresultFile=), analyzeIUPred2A(pathToIUPred2AresultFile=)

aSwitchList <- analyzeSwitchConsequences(
    aSwitchList,
    consequencesToAnalyze = c('intron_retention', 'ORF_seq_similarity', 'NMD_status',
                              'coding_potential', 'domains_identified'),   # + one type per extra annotator imported
    dIFcutoff = 0.1
)
```

| External tool | Purpose | Consequence types it enables | Status |
|----------------|---------|------------------------------|--------|
| CPC2 (standalone, Python 2 code) | Coding vs non-coding | `coding_potential` | Ran on `isoformSwitchAnalyzeR_isoform_nt.fasta`: `python CPC2.py -i <nt.fasta> -o cpc2_result` (needs a Python 2.7 env and its bundled libsvm built) |
| Pfam: `hmmscan --cut_ga --domtblout` against Pfam-A, then `scripts/hmmscan_to_pfamscan.py` | Protein domains | `domains_identified` | Ran (HMMER 3.4 output; below) |
| SignalP 6.0 | Signal peptides | `signal_peptide_identified` | Licence-gated, **not run**; import per `?analyzeSignalP` |
| IUPred2A or NetSurfP | Intrinsically disordered regions | `IDR_identified`; `IDR_type` needs IUPred2A | Licence/registration-gated, **not run**; 2.6.0 source uses `analyzeNetSurfP2()`, while 2.12.0 source uses `analyzeNetSurfP3()`; consult installed help |
| DeepTMHMM | Transmembrane topology | `isoform_topology` (via `analyzeDeepTMHMM()`), not the IDR types | Cloud service, **not run**; read from source only |

Behaviours checked by running:
- **Pfam import.** `analyzePFAM()` rejects a raw `hmmscan --domtblout` file ("more columns than column names") and accepts only a pfam_scan-style table with a clan accession (CL...) in column 15. Convert with the shipped script, which adds clans from Pfam-A.clans.tsv:

  ```bash
  hmmscan --cut_ga --domtblout pfam_domtbl.txt Pfam-A.hmm sequences/isoformSwitchAnalyzeR_isoform_AA.fasta > /dev/null
  curl -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
  python scripts/hmmscan_to_pfamscan.py pfam_domtbl.txt pfam_scanfmt.txt Pfam-A.clans.tsv.gz
  ```
  On the synthetic set (10 genes carrying a planted ubiquitin domain in exon E3): PF00240 found in every isoform containing it and never in the exon-skipped one; 60 domain rows imported; "Domain loss" called 10/10 in the right direction.
- **Missing annotators are an error, not silent.** Requesting a type without its annotation stops: "To test differences in signal peptides, the result of the SignalP analysis must be advailable (sic)"; likewise for CPC2/CPAT (`coding_potential`), Pfam (`domains_identified`) and the installed NetSurfP/IUPred importer (`IDR_*`). Build `consequencesToAnalyze` only from what you imported. A missing file gives "The file(s) ... does not exist".
- **`removeNoncodinORFs` has no default.** `TRUE` drops the ORFs of CPC2-"noncoding" isoforms, which removes PTC-bearing isoforms from the NMD call: on the planted set 8/10 poison isoforms were called NMD-sensitive with `TRUE` (2 were CPC2 noncoding) and 10/10 with `FALSE` or no CPC2. Use `FALSE` when NMD is the question.

The external tools run *outside* R; IsoformSwitchAnalyzeR writes the FASTA files and re-imports the parsed results.
