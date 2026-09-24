---
name: bio-blast-searches
description: Run remote BLAST searches against NCBI servers using Biopython Bio.Blast.NCBIWWW. Use when identifying unknown sequences, finding homologs, picking the correct BLAST program (blastn/blastp/blastx/tblastn/tblastx/psiblast/megablast/dc-megablast), interpreting Karlin-Altschul E-values, avoiding the max_target_seqs trap (Shah 2019), choosing composition-based statistics, or limiting searches by organism. Covers RID lifecycle, database choice (nt/nr/refseq_select/swissprot), word-size and CBS taxonomy.
tool_type: python
primary_tool: Bio.Blast.NCBIWWW
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, NCBI BLAST+ 2.15+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(Bio.Blast.NCBIWWW.qblast)` to check signatures
- CLI: `blastn -version` then `blastn -help`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# BLAST Searches (Remote)

**"Find similar sequences in NCBI's database"** -> Submit a query to NCBI's remote BLAST servers; receive a Request ID (RID); poll for completion; parse the XML hit table. Best for one-off identification of a few sequences at a time (NCBI's queue tolerates roughly one search at a time, a few per minute). For >50 sequences, switch to `local-blast`; for >1000, switch to DIAMOND/MMseqs2 in `remote-homology`.

The two most consequential decisions: **which program** (defines query+target molecule types and word-size defaults) and **which database** (defines the search space and therefore E-value baselines). The third most important: do NOT misuse `max_target_seqs` -- it is an early-termination heuristic, not a "give me the top N hits" filter (Shah et al. 2019).

- Python: `NCBIWWW.qblast(program, db, sequence)` + `NCBIXML.read(handle)` (BioPython)
- CLI: `blastn -remote -db nt -query seq.fa -out hits.xml -outfmt 5` (BLAST+)
- Web: https://blast.ncbi.nlm.nih.gov/Blast.cgi (RID lookup)

## Required Setup

```python
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
```

No API key needed for remote BLAST itself, but NCBI's general rate-limit ethic still applies -- one search at a time, polite waiting, no parallelism.

Always include a FASTA defline in `sequence` (`>id\n...`), even for a throwaway search -- see Failure Modes: "Empty FASTA defline submitted" for the real cost of skipping it.

## Program decision (query vs database molecule)

| Program | Query | Target | Word size default | Use case |
|---|---|---|---|---|
| `blastn` | DNA | DNA | 11 | General DNA similarity |
| `megablast` | DNA | DNA | 28 | High-identity DNA (>=95%) -- PCR primer hits, contamination |
| `dc-megablast` | DNA | DNA | 11 (discontiguous) | Cross-species mRNA (sensitive, gapped) |
| `blastp` | Protein | Protein | 3 (6 also valid) | General protein homology |
| `blastx` | DNA | Protein | 3 | Translated DNA query vs protein DB; ORF discovery |
| `tblastn` | Protein | DNA | 3 | Protein query vs translated DB; find unannotated CDS |
| `tblastx` | DNA | DNA | 3 (both translated) | Most expensive; deep cross-species coding similarity |
| `psiblast` | Protein | Protein | 3 | Iterative PSSM-based remote homology -- see `remote-homology` |

**The misuse to avoid:** using default `blastn` (word=11) for cross-species DNA where `dc-megablast` is the right tool. Or using `megablast` (word=28) for cross-species homology where it will miss every divergent hit. The most-misused BLAST parameter according to literature.

**`qblast()` API note:** `megablast` and `dc-megablast` in the table above are program *concepts*, not literal `program=` strings. `NCBIWWW.qblast(program='megablast', ...)` raises `ValueError: Program specified is megablast. Expected one of blastn, blastp, blastx, tblastn, tblastx` immediately -- confirmed live. Both are requested as `program='blastn'` plus a flag; see "Requesting megablast / dc-megablast" in `references/dna-patterns.md`.

## Database decision (search space)

| Database (`db=`) | Content | Size (2026 approx) | Stable for reproducibility? |
|---|---|---|---|
| `nt` | Non-redundant nucleotide (all GenBank+EMBL+DDBJ) | ~250 GB | NO -- changes daily |
| `nr` | Non-redundant protein | ~300 GB | NO -- changes daily |
| `refseq_select` | One curated rep per species (RNA + protein) | small | YES -- versioned releases |
| `refseq_rna` | RefSeq mRNA | ~10 GB | YES |
| `refseq_protein` | RefSeq protein | small | YES |
| `swissprot` | UniProt Swiss-Prot (reviewed) | small | YES -- monthly releases |
| `pdb` | Protein structures | small | YES |
| `refseq_genomic` | RefSeq genomic | huge | YES |
| `env_nr` / `env_nt` | Environmental (metagenomic) | huge | YES |

**For publication reproducibility, never search `nt` or `nr`** without recording the snapshot date and ideally archiving a frozen copy. Default to `refseq_select` for any cross-species homology question; switch to `nt`/`nr` only when curated coverage is insufficient.

## E-value interpretation (Karlin-Altschul)

E-values scale with database size (Karlin & Altschul 1990 PNAS 87:2264) -- the same alignment against a 100x larger database has a 100x larger E-value. **Cross-database E-value comparison is meaningless; bit-score is database-size normalized and is the right cross-database metric.** For protein remote homology where E is marginal (10^-3 to 10^-1), reach for profile methods: PSI-BLAST, jackhmmer, HHblits, or Foldseek -- see `remote-homology` skill. Full formula, the E-value/bit-score interpretation table, and the homology "twilight zone" threshold: `references/statistics.md`.

## Composition-Based Statistics (CBS)

Compositional bias inflates significance for low-complexity proteins. Default `composition_based_statistics=2` (Yu&Altschul 2005) is correct for most cases; switch to `composition_based_statistics=3` for protein queries under 30 aa, where mode 2 over-corrects (used in the Short peptide search pattern in `references/protein-patterns.md`). For known compositional bias (coiled-coil regions, signal peptides), CBS=2 is appropriate but consider hard-masking with SEG (`filter='S'` in `qblast()`); extreme bias still inflates scores and shows up as many "significant" hits to unrelated low-complexity proteins. Full mode table and mechanism (Yu et al. 2006): `references/statistics.md`.

## The `max_target_seqs` trap

**The misuse**: `max_target_seqs=10` is interpreted as "return the 10 most significant hits". It is not. The flag is an **early termination** parameter that affects which hits the search ever considers, not which it ultimately reports (Shah N, Nute MG, Warnow T, Pop M. (2019) Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows. *Bioinformatics* 35:1613-1614).

**Consequences:**
- Setting `max_target_seqs=10` can return entirely different hits than `max_target_seqs=500` then filtering to top 10 by E-value.
- The "top 10" by E-value as reported may not be the actual top 10.

**Correct pattern:** set `hitlist_size` (Bio.Blast parameter name) large (1000+), then post-filter to the top N by E-value or bit-score in Python.

## Word size, gap costs, and matrix

| Search | Word size | Matrix (protein) | Gap (open, extend) |
|---|---|---|---|
| megablast (high identity DNA) | 28 | n/a | 0, 0 (linear) |
| blastn (sensitive DNA) | 11 | n/a | 5, 2 |
| blastp default | 3 | BLOSUM62 | 11, 1 |
| blastp distant | 2 | BLOSUM45 | 14, 2 |
| Short peptides (<30 aa) | 2 | PAM30 or BLOSUM45 | 9, 1 |

For very short query proteins (e.g. proteomics-identified peptides), BLOSUM45 + word=2 + PAM30 substitution matrix is more sensitive than the default. Use `matrix='PAM30'` for searches against `swissprot`.

## RID lifecycle

| Phase | Server state | Client action |
|---|---|---|
| Submit | RID created, queued | NCBIWWW.qblast() returns handle |
| Running | Queue + compute | Poll status |
| Done | RID + results retained | Fetch XML |
| Expired | RID purged | 24-36h after completion |

`NCBIWWW.qblast()` handles polling internally with a fixed retry interval. For long-running searches (>5 min) or batches, submit and capture the RID, then poll independently to avoid blocking (`scripts/blast_rid.py`, see "Programmatic RID polling" in `references/results-and-rid.md`). The RID is visible at `https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=...` for 24-36 hours.

**Latency is not fixed.** A small blastn against `refseq_select_rna` takes about a minute (62 s measured), but permissive-cutoff, large-`hitlist_size` or low-word-size searches (PAM30 + `word_size=2` + `expect=1000` measured 181 s once and 1141 s (19 min) another time, on an unchanged query pattern) can run many times longer, and any search waits behind a busy NCBI queue (the same small blastn measured 62 s on a quiet day and 33 min on a busy one). A multi-minute wait on those combinations is not a hang: keep polling the RID rather than resubmitting, and only give up on `Status=FAILED`/`UNKNOWN` or after your own timeout.

## Reference Files

Code patterns live in `references/`; read the one that matches the call you are writing.

| File | Read when |
|---|---|
| `references/dna-patterns.md` | Writing a `blastn` call, or requesting megablast / dc-megablast (`megablast=True`, `template_type`) |
| `references/protein-patterns.md` | `blastp` with `entrez_query` organism restriction, or a peptide under 30 aa (PAM30, `gapcosts`) |
| `references/results-and-rid.md` | Saving/re-parsing XML, filtering hits by identity + coverage, or polling a long job by RID |
| `references/statistics.md` | E-value derivation, CBS mode table, twilight zone |

## Failure modes

`max_target_seqs` misinterpretation and cross-database E-value comparison are covered in "The `max_target_seqs` trap" and "E-value interpretation" above.

### Megablast for cross-species
- **Trigger:** Default `megablast` (word=28) on a cross-species DNA query.
- **Mechanism:** Word size 28 requires 28-nt exact match to seed; cross-species mRNA has too much divergence.
- **Symptom:** Reduced or missing cross-species hits, worse for more diverged sequences -- not necessarily zero. Confirmed live: on a real human/mouse/rat query, megablast still found the human and mouse hits (2 of 3 species) and only missed the most-diverged one (rat).
- **Fix:** Use `dc-megablast` (discontiguous) or `blastn` with word=11.

### Reproducibility loss against `nt`/`nr`
- **Trigger:** Manuscript says "BLASTed against nt"; reviewer re-runs 3 weeks later.
- **Mechanism:** Databases change daily; new genomes deposited.
- **Symptom:** Different hit set, different paper conclusions.
- **Fix:** Use `refseq_select` for reproducibility, or record snapshot date + archive subset.

### Server timeout on large queries
- **Trigger:** Multi-megabase query or batch submission.
- **Mechanism:** Remote BLAST has a per-query compute budget.
- **Symptom:** Job stuck in queue, eventually fails.
- **Fix:** Split into smaller queries; or switch to `local-blast` / DIAMOND / MMseqs2.

### Empty FASTA defline submitted
- **Trigger:** Sending `sequence` as a raw string without `>id\n`.
- **Mechanism:** BLAST accepts the anonymous query, but the search runs markedly slower server-side.
- **Symptom:** Two confirmed, reproducible effects, not the one you'd expect from "just a label": (1) `record.query` comes back as the literal placeholder string `'No definition line'`, not `None`/empty; (2) the search itself takes **~12.7x longer** -- confirmed live, 781s vs. 62s for the identical query with vs. without a defline, with Biopython itself raising `BiopythonWarning: BLAST request ... is taking longer than 10 minutes`. Results (alignments, top hit) are otherwise correct.
- **Fix:** Always pass FASTA with a defline (`>id\n...`), or a `SeqRecord`, even for a disposable one-off search -- this is a latency cost, not a cosmetic one, so it is not optional in practice.

## Common errors

| Error / symptom | Cause | Solution |
|---|---|---|
| Stuck > 5 min | Large query or busy queue, or a missing FASTA defline (see Failure Modes: Empty FASTA defline -- confirmed ~12.7x slower) | Submit RID, poll separately; add a defline; or use local |
| URLError / timeout | Network or NCBI maintenance | Retry with backoff; status at status.ncbi.nlm.nih.gov |
| No hits | Wrong program / database type | Verify query and DB molecule types match |
| Empty XML | RID expired | Re-submit; RIDs purge after 24-36h |
| 1000s of low-complexity hits | CBS disabled or extreme bias | CBS=2; consider SEG filter |
| `ValueError: ... Gap existence and extension values of 11 and 1 not supported for PAM30` | Non-BLOSUM62 matrix (e.g. PAM30) without setting `gapcosts` | Set `gapcosts='9 1'` for PAM30 (see word-size/gap-cost table) |
| `ValueError: Program specified is megablast. Expected one of blastn, blastp, blastx, tblastn, tblastx` | `program=` set to `'megablast'`/`'dc-megablast'` directly | Use `program='blastn', megablast=True` (+ `template_type`/`template_length` for dc-megablast) |

## References

- Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. (1990) Basic local alignment search tool. *J Mol Biol* 215:403-410.
- Karlin S, Altschul SF. (1990) Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. *Proc Natl Acad Sci USA* 87:2264-2268.
- Altschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. (1997) Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. *Nucleic Acids Res* 25:3389-3402.
- Yu YK, Gertz EM, Agarwala R, Schaffer AA, Altschul SF. (2006) Retrieval accuracy, statistical significance and compositional similarity in protein sequence database searches. *Nucleic Acids Res* 34:5966-5973.
- Shah N, Nute MG, Warnow T, Pop M. (2019) Misunderstood parameter of NCBI BLAST impacts the correctness of bioinformatics workflows. *Bioinformatics* 35:1613-1614.
- Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. (2009) BLAST+: architecture and applications. *BMC Bioinformatics* 10:421.
- Rost B. (1999) Twilight zone of protein sequence alignments. *Protein Eng* 12:85-94.

## Related Skills

- local-blast - Faster, unlimited local BLAST+ pipelines and database build
- remote-homology - PSI-BLAST, jackhmmer, HHblits, MMseqs2, DIAMOND, Foldseek for distant homology
- ortholog-inference - Reciprocal best hit, OrthoFinder, OMA for orthology
- sequence-io/read-sequences - Load query sequences from FASTA
- entrez-fetch - Fetch full records for BLAST hits
