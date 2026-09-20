'''Build a structural multiple sequence alignment with Foldmason.

Produces both an amino-acid MSA (result_aa.fa) and a 3Di structural MSA (result_3di.fa).
Use --refine-iters > 0 for iterative refinement; --report-mode 1 for an HTML LDDT report,
--report-mode 2 for a machine-readable JSON LDDT report.
Refinement is random unless --refine-seed is set: unseeded --refine-iters 100 runs give different
MSAs (--refine-iters 0 is deterministic), so runs are only comparable with the same seed.
Foldmason writes one MSA row per chain, named <file>_<chain> (a single-chain file keeps its bare
stem), so the row count exceeds the file count for multi-chain files.
'''
# Reference: foldmason 1+ (checked on 4.dd3c235) | Verify CLI flags if version differs

import json
import subprocess
from pathlib import Path

def foldmason_msa(structures, result_prefix, tmp_dir='tmp/', refine_iters=100, report_mode=1, refine_seed=42):
    structure_paths = [str(p) for p in structures]
    cmd = [
        'foldmason', 'easy-msa',
        *structure_paths,
        result_prefix,
        tmp_dir,
        '--refine-iters', str(refine_iters),
        '--refine-seed', str(refine_seed),
        '--report-mode', str(report_mode),
    ]
    subprocess.run(cmd, check=True)
    if report_mode == 1:
        report_path = f'{result_prefix}.html'
    elif report_mode == 2:
        report_path = f'{result_prefix}.json'
    else:
        report_path = None
    return {
        'amino_msa': f'{result_prefix}_aa.fa',
        'structural_msa': f'{result_prefix}_3di.fa',
        'guide_tree': f'{result_prefix}.nw',
        'report': report_path,
    }

def summarize(amino_msa_path):
    '''Return (rows, columns); rows are chains, not input files.'''
    n_seqs = 0
    length = 0
    with open(amino_msa_path) as f:
        for line in f:
            if line.startswith('>'):
                n_seqs += 1
            elif n_seqs == 1:
                length += len(line.strip())
    return n_seqs, length

def per_column_lddt(json_path, n_columns):
    '''Per-column LDDT from a --report-mode 2 JSON: the values are under 'scores' (-1 = not computed).'''
    with open(json_path) as f:
        scores = json.load(f)['scores']
    assert len(scores) == n_columns, f'{len(scores)} scores for {n_columns} MSA columns'
    return scores

if __name__ == '__main__':
    structures = sorted(Path('structures').glob('*.pdb'))
    print(f'Aligning {len(structures)} structures with Foldmason...')

    outputs = foldmason_msa(structures, 'family_msa', refine_iters=100, report_mode=2)
    n_seqs, length = summarize(outputs['amino_msa'])
    print(f'\nMSA: {len(structures)} files -> {n_seqs} chain rows, {length} columns')
    print(f'Amino-acid MSA: {outputs["amino_msa"]}')
    print(f'3Di MSA: {outputs["structural_msa"]}')
    print(f'Guide tree (Newick): {outputs["guide_tree"]}')
    if outputs['report']:
        print(f'LDDT report: {outputs["report"]}')
        scores = per_column_lddt(outputs['report'], length)
        valid = [x for x in scores if x >= 0]
        print(f'Per-column LDDT: {len(valid)} of {len(scores)} columns scored, mean {sum(valid) / len(valid):.3f}')
