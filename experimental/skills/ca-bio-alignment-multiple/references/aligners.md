# Aligners



### Algorithm Selection

MAFFT offers multiple algorithms with explicit accuracy/speed tradeoffs. Selecting the right mode is critical; the difference between L-INS-i and FFT-NS-1 can be the difference between a correct and incorrect downstream phylogeny.

| Algorithm | Flag | Strategy | Best For |
|-----------|------|----------|----------|
| FFT-NS-1 | `--retree 1` | Progressive only | Quick look, >10,000 seqs |
| FFT-NS-2 | `--retree 2` | Progressive + guide tree rebuild | Default balance, 200-10,000 seqs |
| FFT-NS-i | `--maxiterate 1000` | Iterative refinement | Moderate improvement, 200-2,000 seqs |
| G-INS-i | `--globalpair --maxiterate 1000` | Global pairwise + iterative | Sequences alignable over full length, <200 |
| L-INS-i | `--localpair --maxiterate 1000` | Local pairwise + iterative | Single alignable domain amid divergent flanks, <200 |
| E-INS-i | `--genafpair --maxiterate 1000` | Local with generalized affine gaps | Multiple conserved motifs separated by unalignable regions, <200 |
| Auto | `--auto` | Auto-selects based on dataset size | When unsure |

**Decision guide**: If sequences share a single conserved domain (most common case), use L-INS-i. If sequences are globally similar (e.g., ortholog set of similar length), use G-INS-i. If sequences have multiple conserved blocks separated by highly variable linker regions (e.g., multi-domain proteins with variable interdomain regions), use E-INS-i.

**G-INS-i failure mode:** When sequences have long divergent N/C-terminal extensions (signal peptides, intrinsically-disordered regions, isoform-specific tails), the global pairwise distance G-INS-i uses for guide-tree construction is dominated by the extension's mismatch content. The guide tree then mis-clusters sequences by extension similarity rather than core homology. Symptom: the resulting alignment has core conserved domains poorly aligned despite high-identity flanks. Switch to L-INS-i (local pairwise; tolerant of divergent flanks) or pre-process to remove signal peptides/disordered regions before alignment.

### What `--auto` Picks (and Why to Specify Explicitly)

`mafft --auto` silently downgrades the algorithm based on dataset size **and** alignment length (column count) -- verified against the size-check block (`if [ $auto -eq 1 ]`) in the `mafft` wrapper script itself (MAFFT v7.526, `usr/bin/mafft`). Rows are evaluated top-to-bottom; the first match wins.

| Sequences | Columns | `--auto` selects | Equivalent flags |
|-----------|---------|-------------------|-------------------|
| < 100 | < 3,000 | L-INS-i | `--localpair --maxiterate 1000` |
| < 200 | < 1,000 | L-INS-i, but capped to 2 iterations | `--localpair --maxiterate 2` |
| < 500 | < 10,000 | FFT-NS-i | `--retree 2 --maxiterate 2` |
| < 20,000 | any | FFT-NS-2 | `--retree 2 --maxiterate 0` |
| < 100,000 | any | FFT-NS-2, memory-saving guide tree | `--retree 2 --maxiterate 0 --memsavetree` |
| < 200,000 | any | FFT-NS-1, memory-saving guide tree | `--retree 1 --maxiterate 0 --memsavetree` |
| >= 200,000 | < 3,000 | PartTree, local-alignment distance | `--dpparttree` |
| >= 200,000 | >= 3,000 | PartTree, k-mer distance | `--parttree` |

Both sequence count AND column count gate the first three rows -- e.g. 150 sequences of a 1.5 kb gene (150 seqs x ~1,500 columns) fail the `< 100`/`< 3,000` row (too many sequences) and the `< 200`/`< 1,000` row (too many columns), so `--auto` falls through to FFT-NS-i, not L-INS-i, despite being under the old "< 200" threshold. Verified empirically on this machine (MAFFT v7.526): 90 seqs x 200bp -> L-INS-i; 150 seqs x 200bp -> L-INS-i, 2 iterations; 150 seqs x 1,500bp and 120 seqs x 1,100bp -> FFT-NS-i, 2 iterations; 600 seqs x 200bp -> FFT-NS-2. The 100,000 / 200,000 / PartTree rows were confirmed by reading the script only -- running them here wasn't practical.

Separately, the wrapper caps `--maxiterate` at 16 by default no matter what value is requested (`iteratelimit=16` under the default `BAATARI2` parallelization strategy, `usr/bin/mafft` ~line 1512) -- confirmed by running `mafft --localpair --maxiterate 1000` explicitly and seeing "Iterative refinement method (<16)" in the reported strategy, identical to `--auto`'s output for the same input. So hand-specifying `--maxiterate 1000` does not actually raise the iteration count above 16 unless a `--bestfirst`/`--baatari0` parallelization strategy is also set (raises the cap to 254).

For publication-quality phylogenetics, specify the algorithm explicitly so reproducibility audits do not rely on internal threshold heuristics.

### Basic Usage

**Goal:** Run MAFFT on a FASTA file with appropriate algorithm selection.

**Approach:** Invoke MAFFT via command line or subprocess, selecting the algorithm based on dataset characteristics.

```bash
# Highest accuracy for <200 sequences (local pairwise iterative)
mafft --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Good balance for medium datasets
mafft --retree 2 input.fasta > aligned.fasta

# Auto-select algorithm based on dataset size
mafft --auto input.fasta > aligned.fasta

# Protein alignment with specific matrix (default BLOSUM62)
mafft --amino --localpair --maxiterate 1000 input.fasta > aligned.fasta

# DNA alignment (auto-detected, but can be explicit)
mafft --nuc --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Adjust gap penalties (op=gap open, ep=gap extension)
mafft --op 1.53 --ep 0.123 --localpair --maxiterate 1000 input.fasta > aligned.fasta

# Multithreaded
mafft --thread 8 --localpair --maxiterate 1000 input.fasta > aligned.fasta
```

```python
import subprocess

def run_mafft(input_fasta, output_fasta, algorithm='linsi', threads=4):
    algo_flags = {
        'linsi': ['--localpair', '--maxiterate', '1000'],
        'ginsi': ['--globalpair', '--maxiterate', '1000'],
        'einsi': ['--genafpair', '--maxiterate', '1000'],
        'fftns2': ['--retree', '2'],
        'auto': ['--auto'],
    }
    cmd = ['mafft', '--thread', str(threads)] + algo_flags[algorithm] + [input_fasta]
    with open(output_fasta, 'w') as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f'MAFFT failed (exit {result.returncode}):\n{result.stderr}')

run_mafft('sequences.fasta', 'aligned.fasta', algorithm='linsi')
```

MAFFT writes its progress log and errors to stderr, not stdout. Capturing stderr and surfacing the failure message is essential when MAFFT exits non-zero (e.g. encountering ambiguous characters, oversized input, missing libraries) -- otherwise `check=True` raises a CalledProcessError without showing the actionable message.



Pick `-align` (PPP) for peak-accuracy runs up to ~1000 sequences, or `-super5` (mBed clustering + chunked alignment) for thousands to millions. `-super5` is not a "lower quality" mode; both share the same HMM-perturbation single-replicate machinery (`-perm`/`-perturb`), but only `-align` accepts the `-stratified`/`-diversified` ensemble flags directly -- `-super5` rejects them (checked on MUSCLE 5.3.win64: `-super5 ... -stratified` fails with "-stratified not supported").

| Command | Algorithm | Designed for | Output |
|---------|-----------|--------------|--------|
| `-align` (PPP) | Posterior probability progressive (HMM) | <= ~1000 seqs, peak accuracy | Single MSA, or `.efa` ensemble via `-stratified`/`-diversified` |
| `-super5` | mBed-clustering + chunked alignment | Thousands to millions of seqs | Single MSA only per run; no direct `.efa` (combine several `-perm`/`-perturb` runs with `-fa2efa` for an ensemble) |

### Basic Usage

```bash
# Peak accuracy PPP algorithm, <1000 sequences
muscle -align input.fasta -output aligned.fasta -threads 8

# super5 divide-and-conquer for thousands of sequences
muscle -super5 input.fasta -output aligned.fasta -threads 8
```



Use ClustalOmega when datasets reach hundreds of thousands of sequences (the mBed guide tree scales O(N log N)) or for HMM-profile-driven alignment. Prefer MAFFT or MUSCLE5 below that scale.

```bash
# Basic alignment
clustalo -i input.fasta -o aligned.fasta --auto

# Force overwrite output
clustalo -i input.fasta -o aligned.fasta --force

# Specify output format
clustalo -i input.fasta -o aligned.phy --outfmt=phylip

# Use more iterations for better accuracy
clustalo -i input.fasta -o aligned.fasta --iter=5

# Profile-profile alignment (align two existing MSAs)
clustalo --p1 profile1.fasta --p2 profile2.fasta -o merged.fasta

# Add sequences to existing alignment
clustalo -i new_seqs.fasta --profile1 existing.fasta -o updated.fasta

# Multithreaded
clustalo -i input.fasta -o aligned.fasta --threads=8
```
