# ELink code patterns

Loaded from SKILL.md when writing ELink code. Every snippet assumes `from Bio import Entrez`, `import time` and `Entrez.email` set as in SKILL.md "Required Setup".

## Code patterns

### Single source -> single target

**Goal:** Get RefSeq proteins for a single gene.

**Approach:** ELink with explicit `linkname` to restrict to curated set.

**Reference (BioPython 1.83+):**
```python
def gene_to_refseq_proteins(gene_id):
    h = Entrez.elink(dbfrom='gene', db='protein', id=gene_id, linkname='gene_protein_refseq')
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [link['Id'] for link in r[0]['LinkSetDb'][0]['Link']]

print(gene_to_refseq_proteins('672'))  # BRCA1
```

### Batch source -> target (small batch)

**Goal:** Get linked proteins for a list of <200 gene IDs in one call.

**Approach:** Pass the IDs as a list (Biopython sends one `id=` per element); one linkset per input in the response. A comma-joined string would return a single merged linkset (the union), losing which gene each protein belongs to.

**Reference (BioPython 1.83+):**
```python
def batch_gene_protein(gene_ids):
    h = Entrez.elink(dbfrom='gene', db='protein', id=gene_ids, linkname='gene_protein_refseq')
    r = Entrez.read(h); h.close()
    out = {}
    for linkset in r:
        src = linkset['IdList'][0]
        out[src] = [link['Id'] for link in linkset['LinkSetDb'][0]['Link']] if linkset['LinkSetDb'] else []
    return out
```

### Large batch via history server

**Goal:** Link 5,000 gene IDs to proteins without hitting URL-length limits.

**Approach:** EPost the IDs first (chunked), then ELink with `cmd='neighbor_history'` referencing the WebEnv. Downstream EFetch picks up linked IDs from the history server.

**Runnable code:** `examples/chain_links.py`, `link_batch_via_history(dbfrom, db, source_ids, linkname=None, chunk=200)` returns `(WebEnv, QueryKey)`; it unions the per-chunk QueryKeys (see SKILL.md, "Chunked EPost links only the last chunk"). Downstream: `Entrez.efetch(db='protein', WebEnv=we, query_key=qk, retstart=..., retmax=500)`.

### Discover all available links

**Goal:** Before writing a pipeline, enumerate what link tables NCBI exposes for a (dbfrom, source-id) pair.

**Approach:** `cmd='acheck'` returns the full LinkInfo list per source.

**Runnable code:** `examples/discover_links.py`, `discover_links(dbfrom, source_id)` returns `(LinkName, DbTo, MenuTag)` tuples and asserts the LinkInfo schema.

### Chain links (gene -> protein -> structure)

```python
def gene_to_structures(gene_id):
    h = Entrez.elink(dbfrom='gene', db='protein', id=gene_id, linkname='gene_protein_refseq')
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    prot_ids = [l['Id'] for l in r[0]['LinkSetDb'][0]['Link'][:10]]
    time.sleep(0.1 if Entrez.api_key else 0.34)
    h = Entrez.elink(dbfrom='protein', db='structure', id=','.join(prot_ids))
    r = Entrez.read(h); h.close()
    out = []
    for ls in r:
        if ls['LinkSetDb']:
            out.extend(l['Id'] for l in ls['LinkSetDb'][0]['Link'])
    return out
```

### Get neighbor_score for related PubMed articles

```python
def related_pubmed(pmid, top=10):
    h = Entrez.elink(dbfrom='pubmed', db='pubmed', id=pmid,
                     linkname='pubmed_pubmed', cmd='neighbor_score')
    r = Entrez.read(h); h.close()
    if not r[0]['LinkSetDb']:
        return []
    return [(l['Id'], int(l['Score'])) for l in r[0]['LinkSetDb'][0]['Link'][:top]]
```

`Score` is a large raw integer (NCBI's internal relevance magnitude, often in the tens of millions, e.g. `29748057`), not a normalized 0-100 value — never present it as a percentage or small rank.

### BioProject -> SRA runs

For SRA discovery, `pysradb.SRAweb().sra_metadata(prjna, detailed=True)` (see `sra-data`) is the higher-fidelity path — returns SRR accessions directly with run-level metadata in one call. Use ELink only when staying inside Bio.Entrez:

```python
def bioproject_to_sra(prjna):
    # Convert PRJNA to UID first
    h = Entrez.esearch(db='bioproject', term=f'{prjna}[BioProject]')
    r = Entrez.read(h); h.close()
    if not r['IdList']:
        return []
    bp_uid = r['IdList'][0]
    time.sleep(0.1 if Entrez.api_key else 0.34)
    # Link to SRA
    h = Entrez.elink(dbfrom='bioproject', db='sra', id=bp_uid)
    r = Entrez.read(h); h.close()
    return [l['Id'] for l in r[0]['LinkSetDb'][0]['Link']] if r[0]['LinkSetDb'] else []
```
