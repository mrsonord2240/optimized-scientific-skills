# Codon Methods



### When Codon Alignment Is Required

Coding sequences destined for selection analysis (dN/dS with PAML/codeml, HyPhy BUSTED/MEME/aBSREL) **must** be aligned respecting codon boundaries. Standard nucleotide MSA tools do not preserve reading frames and produce systematically incorrect dN/dS estimates -- Fletcher & Yang (2010 MBE) showed conventional aligners cause false-positive signals in the branch-site test of positive selection even on clean simulated data.

### Codon-Alignment Tool Decision Tree

| Input cleanliness | Recommended tool | Notes |
|------------------|------------------|-------|
| Clean orthologs (no frameshifts, no internal stops) | MAFFT-protein + PAL2NAL (Suyama, Torrents & Bork 2006 NAR) | Fastest; standard PAML pipeline input |
| Recently duplicated paralogs (indel-rich) | PRANK +F codon (Loytynoja & Goldman 2008 Science) | Phylogeny-aware indel model; fewest false-positive selection calls (Fletcher & Yang 2010) |
| Frameshifts, pseudogenes, error-prone assemblies | MACSE v2 `alignSequences -fs <cost>` (Ranwez et al 2018 MBE) | Frameshift-tolerant; preserves reading frame across sequencing errors |
| Mixed dataset with some bad genes | OMM_MACSE pipeline (Scornavacca, Belkhir, Lopez et al 2019; Ranwez group) | MACSE non-homologous-fragment trimming + MAFFT prealignment + MACSE frameshift refinement + HMMcleaner |
| HyPhy-grade quality (BUSTED, MEME, aBSREL input) | HyPhy `pre-msa.bf` / `post-msa.bf` (Pond lab; Kosakovsky Pond et al) | Strips stop codons, runs MSA at protein level, threads back, validates frames |

### PAL2NAL: Protein-Guided Codon Alignment

**Goal:** Thread a nucleotide coding sequence alignment onto a protein alignment to preserve reading frame.

**Approach:** Align protein sequences first (higher sensitivity), then use PAL2NAL to map the protein alignment back to codons. Suyama et al 2006 NAR; standard codeml input pipeline.

```bash
mafft --localpair --maxiterate 1000 proteins.fasta > proteins_aligned.fasta
pal2nal.pl proteins_aligned.fasta codons.fasta -output fasta > codons_aligned.fasta
pal2nal.pl proteins_aligned.fasta codons.fasta -output paml > codons_aligned.phy
```

**Non-standard genetic codes:** PAL2NAL by default uses the standard code (NCBI table 1). For mitochondrial (vertebrate=2, yeast=3, invertebrate=5), ciliate macronuclear (=6, UAA/UAG=Gln), or other non-standard codes, the protein-to-codon mapping mismatches and PAL2NAL silently produces wrong codon assignments. Specify explicitly with `-codontable N`:

```bash
pal2nal.pl proteins_aligned.fasta codons.fasta -output paml -codontable 2 > codons_mt.phy
```

For datasets mixing genetic codes (e.g. nuclear + mitochondrial CDS in one tree), translate each lineage with its own code BEFORE protein alignment, never apply a single code globally. See NCBI Translation Tables for the full numbering.

### PRANK +F Codon Mode

**Goal:** Align coding sequences under a phylogeny-aware indel model that does not over-collapse insertions.

**Approach:** PRANK +F (Loytynoja & Goldman 2008 Science) treats indels as evolutionary events on a tree, distinguishing insertions from deletions. Slower than MAFFT but recommended for selection-analysis prep.

```bash
prank -d=codons.fasta -o=prank_aligned -codon -F
```

The `+F` flag enforces "fewer false insertions"; without it PRANK behaves more like a conventional aligner. For dN/dS analysis under PAML M2a/M8, PRANK +F is the most conservative input choice.

### MACSE v2: Frameshift-Tolerant Codon Alignment

**Goal:** Align coding sequences directly at the codon level while tolerating frameshifts and internal stops.

**Approach:** MACSE scores based on amino acid translation while operating on DNA. The v2 toolbox (Ranwez et al 2018 MBE) ships a sub-program selector via `-prog`; the eight most-used sub-programs are listed below (run `macse -help` for the full set in the installed release):

| Sub-program | Purpose |
|-------------|---------|
| `alignSequences` | Align coding sequences (`-fs <cost>` sets frameshift penalty) |
| `enrichAlignment` | Add new sequences to existing codon alignment |
| `refineAlignment` | Refine an alignment with MACSE scoring |
| `splitAlignment` | Split into sub-alignments (e.g. by guide tree clades) |
| `exportAlignment` | Convert between MACSE-internal and standard formats |
| `trimAlignment` | Position-aware trimming preserving codon structure |
| `trimNonHomologousFragments` | Detect and remove non-homologous regions per sequence |
| `reportGapsAA2NT` | Map protein-level gap removals back to nucleotide alignment |

```bash
java -jar macse_v2.jar -prog alignSequences -seq coding_seqs.fasta \
    -out_NT aligned_nt.fasta -out_AA aligned_aa.fasta -fs 30

java -jar macse_v2.jar -prog enrichAlignment -align existing.fasta \
    -seq new_seqs.fasta -out_NT updated.fasta
```

### OMM_MACSE: Recommended Pipeline Wrapper

OMM_MACSE (Ranwez group, used in OrthoMaM v10+) chains MACSE non-homologous-fragment trimming, MAFFT pre-alignment, MACSE v2 frameshift-aware refinement, and HMMcleaner -- the recommended pipeline for genome-scale ortholog datasets where some genes will contain frameshifts. See `references/specialist-tools.md` for the pipeline description and invocation command.

### HyPhy pre-msa.bf / post-msa.bf

For HyPhy-grade dN/dS analyses (BUSTED, MEME, aBSREL, RELAX), Pond lab's standard workflow strips stop codons, aligns at the protein level, then threads back and validates frames. See `references/specialist-tools.md` for the three-command sequence.
