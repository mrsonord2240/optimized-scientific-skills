## Tool Taxonomy

| Tool (lang) | BEAST `[&...]` | NHX | phyloXML richness | When |
|-------------|----------------|-----|-------------------|------|
| treeio (R) | YES, structured (`read.beast`/`read.mrbayes`/`read.iqtree`) | YES | via tidytree/ggtree | the default whenever annotations matter; feeds ggtree |
| DendroPy (Py) | YES, structured (`.annotations`) | YES | no | Python work needing metadata, posterior sets, tree distances, conversion with annotations |
| ETE3/ETE4 (Py) | partial (custom features) | YES, native | no | reconciliation, NHX round-trips, programmatic node features |
| Bio.Phylo (Py) | NO key-value parse (opaque `.comment`) | no | YES, richest | general pipelines, conversion among its 5 formats, phyloXML annotation |
| ape (R) | NO (drops `[&...]`) | no | no | fast topology/branch-length analysis; pair with treeio for annotated files |

References: Bio.Phylo Talevich 2012; DendroPy Sukumaran 2010; ETE3 Huerta-Cepas 2016; ape Paradis 2019; treeio Wang 2020.

One-line decision rule: if the file came from BEAST, MrBayes, or RevBayes or carries `[&...]`/`&&NHX` that matters, read it with treeio (R) or DendroPy (Py); otherwise Bio.Phylo (Py) or ape (R) is fine. Never route a BEAST MCC tree through `Bio.Phylo` or `ape::read.nexus` when the HPDs are needed.

## Format Capability (What Each Format Can Hold)

| Capability | Newick | NHX | Nexus (BEAST-annotated) | phyloXML | NeXML |
|------------|--------|-----|-------------------------|----------|-------|
| Topology + branch lengths | yes | yes | yes | yes | yes |
| One support value | ambiguous slot | `B=` tag | comment key | typed `<confidence>` | typed meta |
| Multiple supports per node | no (single slot) | tags | comment keys | yes (n elements) | yes |
| Posterior / HPD intervals | no | no | yes (`_95%_HPD={}`) | via property | via meta |
| Per-branch rates / dates | no | tags | yes | yes | yes |
| Taxonomy (NCBI id/rank) | no | `S=`/`T=` | no | yes, typed | yes |
| Schema-validated | no | no | loose | yes (XSD) | yes (XSD) |

Newick/NHX are compact and grep-able; phyloXML/NeXML are verbose XML (often 5-20x larger) but typed and validatable -- prefer XML for archiving/exchange where machine-checkable semantics matter, Newick/Nexus for pipeline interchange. BEAST/FigTree metadata rides inside a Nexus (or Newick) `[&...]` comment, so "it is a Nexus file" says nothing about whether annotations survive -- only the parser does.

## Reading, Writing, and Converting (Bio.Phylo)

**Goal:** Move trees between formats and inspect them without assuming a single tree or losing annotations.

**Approach:** Use `Phylo.read` for exactly one tree and `Phylo.parse` for many (posterior sets); use `Phylo.convert` only among formats of equal or greater capability; check `.confidence` vs `.name` to confirm support was read into the right slot.

```python
from Bio import Phylo

tree = Phylo.read('tree.nwk', 'newick')          # exactly one tree; raises if 0 or >1
posterior = list(Phylo.parse('run.trees', 'nexus'))   # many trees: posterior/bootstrap set
Phylo.write(tree, 'tree.xml', 'phyloxml')        # phyloXML is Bio.Phylo's richest format

Phylo.convert('tree.nex', 'nexus', 'tree.nwk', 'newick')   # WARNING: [&...] is written as an escaped [\[&...\]] comment that DendroPy/treeio cannot parse -- annotations are effectively corrupted

for clade in tree.get_nonterminals():
    print(clade.confidence, clade.name)          # confirm the support landed in .confidence, not .name
    if clade.confidence is None and clade.name and '/' in clade.name:   # IQ-TREE -B + --alrt: '98.5/100' stays in .name
        sh_alrt, ufboot = (float(v) for v in clade.name.split('/')[-2:])   # last field = UFBoot
```

Supported format strings: `newick`, `nexus`, `phyloxml`, `nexml`, `cdao` (`cdao` raises `KeyError` unless `pip install rdflib`, and does not keep confidences). Colors and branch widths persist only in phyloXML.

Encoding: on Windows, path-based `Phylo.read`/`Phylo.write` (and DendroPy `path=`) use the locale codec (cp1252), garbling non-ASCII tip names. Pass explicit handles for both directions: `with open(p, encoding='utf-8') as fh: Phylo.read(fh, 'newick')`, and the same with `'w'` for writing.

## Preserve BEAST/MrBayes Annotations Before Down-Converting

**Goal:** Keep posteriors, HPD intervals, and rates when a downstream tool wants plain Newick.

**Approach:** Read with an annotation-aware parser, extract the numbers into a side table that travels with the analysis, and only then write a stripped topology -- never down-convert first.

```python
import dendropy

tree = dendropy.Tree.get(path='mcc.tree', schema='nexus', extract_comment_metadata=True)
for node in tree:
    if node.annotations.get_value('posterior') is not None:
        post = float(node.annotations.get_value('posterior'))   # DendroPy returns strings; cast
        hpd = [float(v) for v in node.annotations.get_value('height_95%_HPD')]   # raw BEAST key; treeio typically exposes it as height_0.95_HPD (exact name varies by source program and version -- introspect the columns)
        # persist post/hpd to a side table keyed by the clade before any conversion
tree.write(path='topology.nwk', schema='newick', suppress_annotations=True, suppress_rooting=True)   # intentional, after extraction; suppress_rooting drops the leading [&R]
```

In R the equivalent is treeio `read.beast('mcc.tree')` then `get.data()` / `as_tibble()`, feeding ggtree (tree-visualization); `write.beast()` re-serializes with annotations intact.

## Per-Method Failure Modes

### BEAST/MrBayes MCC to Plain Newick Erases the Credible Intervals
**Trigger:** `Phylo.convert`, `ape::read.nexus` + `write.tree`, or any "just give me the topology" step on an annotated tree.
**Mechanism:** The HPDs, posteriors, and rates live only in the `[&...]` comments; stripping parsers (`ape::read.nexus`) discard them, and Bio.Phylo writes them as escaped `[\[&...\]]` comments that DendroPy and treeio no longer parse.
**Symptom:** The output plots fine but the credible intervals are gone, irrecoverable without re-running the MCMC.
**Fix:** Read with treeio `read.beast` or DendroPy `extract_comment_metadata=True`; extract the numbers to a side table; keep the original `.tree` as the source of truth.

### Support Value Read as a Node Name (or Truncated)
**Trigger:** Parsing a tree whose internal-node slot holds a bootstrap, a posterior, a clade name, or IQ-TREE's `SH-aLRT/UFBoot` dual value.
**Mechanism:** The Newick grammar gives one slot for all of these; the parser must be told which it is. Bio.Phylo keeps a non-numeric label such as `98.5/100` whole in `.name` and leaves `.confidence` None.
**Symptom:** Supports appear as `.name` strings, and any code reading `.confidence` silently gets None.
**Fix:** Know what wrote the file; in Bio.Phylo inspect `.confidence` vs `.name` and split `name.split('/')` (SH-aLRT, UFBoot last) as in the snippet above; in treeio use `read.iqtree`/`read.raxml`, which split dual support correctly.

### Whitespace, Underscore, or Non-ASCII Taxon Names
**Trigger:** Tip names with spaces, parentheses, commas, or accented characters; reliance on the Newick underscore-space convention.
**Mechanism:** Naive CLI tools split unquoted spaces, and underscore-to-space auto-conversion silently desyncs tip labels from a metadata join key.
**Symptom:** Downstream tools error or a metadata merge matches nothing.
**Fix:** Read and write through UTF-8 handles; sanitize to `[A-Za-z0-9_.]` before ape or ete3 ever see the file (ape reads `'O''Brien isolate 7'` as `NA`; ete3 raises on quoted names), keep a label map, and round-trip-test the labels against the metadata table before any join.

### Nexus TRANSLATE-Table Desync and Rooted/Unrooted Confusion
**Trigger:** Hand-editing or merging Nexus tree blocks; assuming topology shape implies rootedness.
**Mechanism:** The integer-to-name TRANSLATE map can decouple from the tree and silently relabel tips; Newick does not flag rootedness (a basal trifurcation conventionally signals unrooted, but tools disagree), while Nexus carries an explicit `[&R]`/`[&U]`.
**Symptom:** Tips are mislabeled after a merge, or a rooting-sensitive analysis runs on the wrong assumption without erroring.
**Fix:** Parse with a translate-aware reader (treeio/DendroPy/ape apply it); verify tip-label sets match across merged trees; set rootedness explicitly rather than trusting topology shape.

## Quantitative and Practical Notes

| Item | Guidance | Why |
|------|----------|-----|
| Support-value scales | bootstrap/UFBoot in [0,100], posterior in [0,1], SH-aLRT in [0,100] | a number is meaningless without knowing which test produced it; preserve provenance, not just the value |
| Multi-tree files | use `Phylo.parse` / DendroPy `TreeList` / treeio `read.beast`; check object length | a single-tree reader on a `.trees` posterior returns only the first or errors |
| Round-trip test | read -> write -> read and diff the annotations, not just the topology | topology almost always survives and gives false confidence |
| Posterior set vs MCC | `.trees` is the full posterior; `.tree`/`.mcc` is the single annotated summary | `read.beast` on a full posterior is huge; usually the MCC is wanted |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| HPD bars missing after conversion | converted a BEAST tree to Newick | extract annotations with treeio/DendroPy first |
| `Phylo.read` raises on a `.trees` file | multiple trees in the file | use `Phylo.parse` and iterate |
| Bootstrap values show up as taxon names | node-label slot read as `.name` | set/inspect confidence parsing; use a software-specific reader |
| Metadata join matches nothing | underscore/space relabeling of tips | sanitize and round-trip-test labels before joining |
| Parser errors on `[` | strict parser chokes on FigTree comment | strip comments only after extracting needed metadata |
| `AssertionError: Two string taxonomies?` from `Phylo.read` on a MrBayes `.con.tre` | Bio.Phylo cannot parse MrBayes' annotated consensus layout | read with DendroPy (`extract_comment_metadata=True`) or treeio `read.mrbayes` |
| NeXML tips named `d7`, `d8`... | Bio.Phylo's NeXML reader names tips by otu id, not label | check tip names; read NeXML from other tools with DendroPy |
| NeXML written by `Phylo.write(..., 'nexml')` has no taxonomy, and drops confidence for phyloXML-sourced trees | the writer never touches `.taxonomies`, and only checks the singular `.confidence` -- phyloXML clades carry support in `.confidences` (a typed list: bootstrap, probability, ...), which stays `None` and gets skipped; a plain single-value `.confidence` (e.g. from Newick) IS written, as a `cdao:has_Support_Value` edge property | write annotated NeXML with DendroPy instead -- it serializes whatever it read (posterior, HPD, rate, taxonomy, ...) as typed `<meta>` elements, not just one hardcoded slot |
| Non-ASCII tip garbled (`CercopithÃ¨que`) | path-based read/write used the Windows locale codec | open files with `encoding='utf-8'` for read and write |
