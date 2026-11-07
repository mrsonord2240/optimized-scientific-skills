## Algorithmic Taxonomy

| Family | Unit | Reference tools | Fails when |
|--------|------|-----------------|------------|
| Event-based | Pre-defined SE/A5SS/A3SS/MXE/RI events from annotation | rMATS-turbo, SUPPA2, VAST-TOOLS | Event isn't in annotation; complex multi-junction events split arbitrarily; AFE/ALE confounded with splicing |
| LSV-based | Local Splice Variations at single source/target nodes | MAJIQ V3 | Memory-constrained environments; cohorts smaller than ~3 reps; non-academic users (license) |
| Junction-cluster | Annotation-free intron clusters by shared splice sites | leafcutter, leafcutter2 | Undersampled clusters lose power; topology biologically uninterpretable for novel events |
| Splice-graph | Graph nodes (non-overlapping exonic regions) | Whippet, Shiba | Whippet maintenance status uncertain since ~2022; complex multi-exon graphs |
| Coverage-aware (IR) | Intron body coverage + flanking junctions | IRFinder, IRFinder-S, S-IRFindeR, iREAD | Confounded by overlapping exons, repeats, low mappability regions |
| Isoform-based | Transcript abundance via EM | Salmon/kallisto + tximport | Salmon EM uncertainty propagates; many similar isoforms (TTN, MAPT) become indistinguishable |

The agent's first decision is **which family** the question requires, not which tool. Switching family is the response to a tool failing within a family — switching tool within a family rarely fixes systematic blind spots.

## Event Taxonomy (Beyond Standard SE/A5SS/A3SS/MXE/RI)

| Class | Code | Biology | Detection caveat |
|-------|------|---------|-------------------|
| Skipped exon (cassette) | SE | Default cassette-exon AS | Most common; well-handled by all tools |
| Alternative 5' splice site | A5SS | Alternative donor; intron 5' end varies | Sign convention tool-specific (see below) |
| Alternative 3' splice site | A3SS | Alternative acceptor; intron 3' end varies | Sensitive to BPS cancer mutations (SF3B1) — cryptic 3'ss ~10-30nt upstream |
| Mutually exclusive exons | MXE | Two exons paired, one included | Tool implementations vary; verify which "form 1" is yours |
| Retained intron | RI | Whole intron retained in mature mRNA | Junction-only quant systematically underdetects; needs IRFinder |
| **Microexon** | (sub-SE) | 3-27 nt; neural-enriched, SRRM4-regulated | Missed by default aligner anchor lengths (>=20-30 nt); needs VAST-TOOLS, MicroExonator, or long-read |
| **Exitron** | n/a | Intronic region within an annotated CDS exon | Mis-classified as A5SS/A3SS by most tools; use ExitronFinder, ScanExitron |
| **Alternative first exon (AFE)** | AFE | Alternative TSS / promoter use | Promoter-driven, NOT spliceosomal; confirm with FANTOM CAGE before reporting as splicing |
| **Alternative last exon (ALE)** | ALE | Alternative cleavage/polyadenylation | APA-driven, NOT spliceosomal; confirm with 3'-end-seq |
| **Detained intron (DI)** | (RI subtype) | Nuclear-retained on mature mRNA, regulated by Clk kinases | Distinct from cytoplasmic NMD-targeted RI (Boutz 2015 *Genes Dev*); requires fractionation to confirm |
| **Recursive splicing** | (intron subtype) | Long introns >50kb spliced via internal ratchet points | Sibley 2015 *Nature*; only detectable with long-read or nascent RNA-seq |

Tool-agnostic taxonomy reference: Wang 2008 *Nature*; Vaquero-Garcia 2016 *eLife* (LSV); Tapial 2017 *Genome Res* (VastDB).

## Tool Selection Matrix

| Tool | Best for | Input | Strengths | Fails when |
|------|----------|-------|-----------|------------|
| rMATS-turbo | Standard SE/A5SS/A3SS/MXE/RI in annotated organism, n>=3 | BAM + GTF | Fast, well-calibrated at n>=3, novel SS support | Junction read imbalance; novel multi-junction events; underdetects RI and microexons |
| SUPPA2 | Quick PSI from existing TPM, pilot analysis | Salmon/kallisto TPM + GTF | No alignment; fastest | Annotation-bound; high FDR (15-30%) at n<=2 vs n<=2 |
| MAJIQ V3 | Complex events, heterogeneous cohorts (HET) | BAM + GFF3 | Bayesian posterior PSI, complete LSV semantics | High memory on large cohorts; academic license; complex LSV interpretation needs care |
| leafcutter | Novel junction discovery, sQTL, low-memory environments | BAM (regtools junctions) | Annotation-free; low memory; suited to unannotated organisms | Sensitive to read depth; cluster topology arbitrary for complex multi-junction events |
| Shiba (2025) | Low-coverage / few-replicate designs; junction imbalance correction | BAM + GTF | Best calibration in own benchmark at n=2 vs n=2 (not executed here) | New (2025); limited community calibration |
| VAST-TOOLS | Cross-species comparative AS, microexons | FASTQ | VastDB orthology, ExOrthist co-tool (not executed here; 6.7 GB VASTDB) | Limited to species in VastDB |
| Whippet | Laptop-scale exploratory | FASTQ | Fast splice-graph-based PSI | Reduced active development since ~2022; underperforms on complex topologies |
| IRFinder / IRFinder-S | Intron retention specifically | FASTQ | Coverage + junction integration; IRFinder-S adds CNN-based artifact filtering (IRFinder 1.3.1 executed; IRFinder-S 2.0 not) | IR-only; not for cassette events |
| S-IRFindeR | Replicate-stable IR ratio | BAM | Stable IR ratio metric | IR-only; less integrated than IRFinder-S |

Methodology evolves; verify benchmarks (Olofsson 2023 *Biochem Biophys Res Commun*; Kubota 2025 *NAR*; Tran 2025 *WIREs RNA*) and tool docs before committing. Default 2026 recommendation: run rMATS-turbo + leafcutter and reconcile; add MAJIQ V3 for complex events / heterogeneous cohorts; switch to Shiba for n=2 vs n=2.

## PSI Definition and Effective Length Normalization

For a cassette exon, naive PSI ignores that the inclusion isoform contains more positions where a junction read can map than the skipping isoform. rMATS reports `IncFormLen` and `SkipFormLen` and computes (SE shown). In the JC files the maxima are 2*(readLength - 1) and readLength - 1 (98 and 49 at `--readLength 50`; 148 and 74 at 75), reached only for exons at least read-length long; short exons shorten them (on real chrX data 719 of 958 SE rows were 148/74), so use the per-event `IncFormLen`/`SkipFormLen` rMATS reports in each row, not the constants; the JCEC files add exon-body positions to `IncFormLen` (149 for the 100 nt planted exon):

```
PSI = (IJC / IncFormLen) / (IJC / IncFormLen + SJC / SkipFormLen)
```

where IJC = inclusion junction counts, SJC = skipping junction counts. **Skipping this normalization biases PSI**: planted counts IJC 80, SJC 10 give naive 0.889 but length-normalized 0.800 (the true value). SUPPA2 derives PSI from transcript TPMs, which the upstream Salmon/kallisto already accounts for. leafcutter operates on intron usage proportions within a cluster (a different statistic).

For long-read data, every read carries full isoform identity — effective-length normalization becomes unnecessary because each read counts as one isoform.

## Sign Conventions for Alternative Splice Sites

| Tool | A5SS interpretation | A3SS interpretation | "Inclusion" direction |
|------|----------------------|----------------------|------------------------|
| rMATS | "long" form = donor downstream of alternative donor | "long" form = acceptor upstream of alternative acceptor | PSI > 0 = more long form |
| SUPPA2 | Same as rMATS (long = inclusion of additional exon body) | Same | PSI > 0 = more long form |
| VAST-TOOLS | Encoded in event ID (`_D1` vs `_D2`) | Encoded in event ID (`_A1` vs `_A2`) | Document the chosen reference |
| MAJIQ | Per-junction within LSV; explicit donor/acceptor naming in VOILA | Same | PSI per junction in the LSV |

Always record which alternative form ΔPSI > 0 corresponds to in publication-grade reporting. Confusion is the most common reviewer comment for AS papers.
