# pysam Consensus, Comparison and Header Dict

Moved from `SKILL.md` (2026-09-21); `build_consensus` and `compare_to_ref` are now `scripts/pysam_consensus.py`. Read for the pysam pileup consensus (teaching only; use `samtools consensus` for real work), the per-position comparison to the reference, and the header dict for writing a BAM.

### Generate Simple Consensus
```python
import pysam
from collections import Counter

def consensus_at_position(bam, chrom, pos):
    bases = Counter()
    for pileup in bam.pileup(chrom, pos, pos + 1, truncate=True):
        if pileup.pos == pos:
            for read in pileup.pileups:
                if not read.is_del and not read.is_refskip:
                    bases[read.alignment.query_sequence[read.query_position]] += 1
    if bases:
        return bases.most_common(1)[0][0]
    return 'N'

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    consensus = consensus_at_position(bam, 'chr1', 1000000)   # 0-based position
    print(f'Consensus at chr1:{1000000 + 1} = {consensus}')
```

### Build Consensus Sequence (Pedagogical Only)

The Python majority-vote consensus below is illustrative, NOT production. `samtools consensus` is Bayesian, quality-aware, and platform-aware; majority vote weights every base equally and produces wrong calls on low-coverage / low-quality regions. Use for teaching pileup iteration mechanics; use `samtools consensus` for any real consensus.

`pileup()` filters before the vote (pysam 0.24.1 defaults): bases with quality < 13, unmapped / secondary / QC-fail / duplicate reads, and orphan reads are dropped, and overlapping mates are counted once. `max_depth` defaults to 8000 (deeper columns are subsampled), so raise it.

`build_consensus` returns exactly `end - start` characters, so index `i` is reference position `start + i`. `pileup()` skips uncovered columns; the function starts from all-`N` and fills only the columns it sees. Building the string by appending per pileup column shifts everything after the first coverage gap (581 false differences vs 1 true on the real chr22 slice). A column deleted in every read counts no base, so it is `N`; insertions are ignored.

Code: `scripts/pysam_consensus.py` (`build_consensus`, importable; window is 0-based half-open):
```bash
python scripts/pysam_consensus.py consensus input.bam chr22 1951 4617 --min-depth 3   # prints the consensus string
```

### Compare Consensus to Reference (Python)
`compare_to_ref` in `scripts/pysam_consensus.py` returns `[(1-based position, ref base, consensus base)]` for called bases that differ from the reference; the CLI prints them tab-separated:
```bash
python scripts/pysam_consensus.py compare input.bam reference.fa chr22 1951 4617 --min-depth 3
```
A reference `N` or IUPAC code never equals a called base, so those positions are listed as differences (a 30-base `N` run gave 30 entries). Ties (50/50 columns) go to the first base counted; `samtools consensus` calls them `N` (or an IUPAC code with `--ambig`) and weights bases by quality, so expect a few different calls at het columns and at shallow, low-quality columns (chr22 slice 1952-4617 at `-d 3`: 1 difference by majority vote and by `-m simple --call-fract 0.5 --min-BQ 13`, 2 by the default Bayesian mode; at the default `-d 1`: 5 and 4).

### Header Dict for Writing a BAM (not a .dict file)
`pysam.AlignmentFile(..., 'wb', header=header)` takes this dict. It has no `M5`, so it is not a sequence dictionary: use `samtools dict` for that.
```python
import pysam

def create_dict_header(fasta_path):
    header = {'HD': {'VN': '1.6', 'SO': 'unsorted'}, 'SQ': []}

    with pysam.FastaFile(fasta_path) as ref:
        for name in ref.references:
            length = ref.get_reference_length(name)
            header['SQ'].append({'SN': name, 'LN': length})

    return header

header = create_dict_header('reference.fa')
for sq in header['SQ'][:5]:
    print(f'{sq["SN"]}: {sq["LN"]:,} bp')
```
