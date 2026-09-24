# PubMed retrieval

Read when fetching PubMed records (abstract, MEDLINE, XML with MeSH/grants/PMC ID).

Assumes the Required Setup block in `SKILL.md` (`Entrez.email`, imports, `expect_start`).

### pubmed

| rettype | retmode | Returns | Use when |
|---|---|---|---|
| `abstract` | `text` | Title + authors + abstract | Reading abstracts |
| `medline` | `text` | MEDLINE flat | Parsing with `Bio.Medline` |
| `xml` | `xml` | Full PubMed XML | Programmatic — get MeSH, grants, PMC link |
| (omitted) | (omitted) | Defaults to XML | EFetch default for pubmed is XML — pass `retmode='xml'` explicitly for clarity |

### Pull PubMed with structured MeSH

**Goal:** Get MeSH terms and grant information that aren't in the abstract format.

**Approach:** `rettype='xml'` and walk the PubmedArticle structure defensively.

**Reference (BioPython 1.83+):**
```python
def pubmed_full(pmid):
    h = Entrez.efetch(db='pubmed', id=pmid, retmode='xml')
    records = Entrez.read(h); h.close()
    article = records['PubmedArticle'][0]
    citation = article['MedlineCitation']
    mesh = [m['DescriptorName'] for m in citation.get('MeshHeadingList', [])]
    title = citation['Article']['ArticleTitle']
    # ArticleIdList entries are StringElement (a str subclass with .attributes), not
    # {'#text': ...} dicts, on Biopython 1.88 -- id['#text'] raises TypeError. Use str(id).
    ids = article.get('PubmedData', {}).get('ArticleIdList', [])
    pmc_id = next((str(i) for i in ids if i.attributes.get('IdType') == 'pmc'), None)
    return {'pmid': pmid, 'title': title, 'mesh': mesh, 'pmc_id': pmc_id}
```
