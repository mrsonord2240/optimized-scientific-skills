# Per-database link catalog

Loaded from SKILL.md when a request needs the linknames for a specific (`dbfrom`, `db`) pair.

## Per-database link catalog (curated subset)

### gene

| Target | Common linknames | Notes |
|---|---|---|
| protein | `gene_protein`, `gene_protein_refseq`, `gene_protein_swissprot` | RefSeq is the safe default |
| nuccore | `gene_nuccore`, `gene_nuccore_refseqrna`, `gene_nuccore_refseqgene` | `refseqrna` for mRNA, `refseqgene` for the curated gene region |
| pubmed | `gene_pubmed`, `gene_pubmed_rif` | RIF is curated and high-quality |
| homologene | `gene_homologene` | Deprecated 2014, data still queryable but new entries stopped -- prefer Ensembl Compara or OrthoFinder for current orthology |
| snp | `gene_snp` | dbSNP entries in gene region |
| clinvar | `gene_clinvar` | Clinical variants |
| omim | `gene_omim` | Disease associations |

### nuccore / nucleotide

| Target | Common linknames |
|---|---|
| protein | `nuccore_protein`, `nuccore_protein_refseq` |
| gene | `nuccore_gene` |
| taxonomy | `nuccore_taxonomy` |
| biosample | `nuccore_biosample` |
| sra | `nuccore_sra` |
| pubmed | `nuccore_pubmed`, `nuccore_pubmed_refseq` |

### protein

| Target | Common linknames |
|---|---|
| nuccore | `protein_nuccore`, `protein_nuccore_cds`, `protein_nuccore_mrna` |
| gene | `protein_gene` |
| structure | `protein_structure` |
| cdd | `protein_cdd` (conserved domains) |
| pubmed | `protein_pubmed` |

### pubmed

| Target | Common linknames |
|---|---|
| pubmed | `pubmed_pubmed`, `pubmed_pubmed_citedin`, `pubmed_pubmed_refs` |
| gene | `pubmed_gene`, `pubmed_gene_rif` |
| protein | `pubmed_protein` |
| nuccore | `pubmed_nuccore` |
| gds | `pubmed_gds` (GEO datasets cited in paper) |
| sra | `pubmed_sra` |

### bioproject

| Target | Common linknames |
|---|---|
| biosample | `bioproject_biosample` |
| sra | `bioproject_sra` |
| pubmed | `bioproject_pubmed` |
