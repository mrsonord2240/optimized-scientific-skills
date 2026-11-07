'''Empirical p-value via shuffling for pairwise alignment.

Reference: BioPython 1.83+ | Verify API if version differs
Reference: Pearson 2013 Curr Protoc Bioinf 3.1; Altschul-Erickson 1985 MBE.

Two shuffle strategies:
- Mononucleotide / mono-residue shuffle: preserves composition only
- Dinucleotide / di-residue shuffle: preserves composition + transition frequencies
  (Altschul-Erickson 1985 algorithm; better null for biased compositions)

For protein, mono-residue shuffle is the standard; dinucleotide is mainly
relevant for DNA where local context is meaningful (preserve='di' below, a pure-Python
Altschul-Erickson shuffle; no ushuffle build needed).
'''
import random
from Bio.Align import PairwiseAligner, substitution_matrices

def dinuc_shuffle(seq):
    '''Uniform random sequence with the same dinucleotide counts, first and last letter
    (Altschul & Erickson 1985: random Eulerian path via a random last-exit tree).'''
    seq = str(seq)
    if len(seq) < 3:
        return seq
    succ = {}
    for a, b in zip(seq, seq[1:]):
        succ.setdefault(a, []).append(b)
    last = seq[-1]
    while True:   # draw a last-exit edge per letter; keep only draws whose edges all lead to `last`
        exit_edge = {v: random.choice(nxt) for v, nxt in succ.items() if v != last}
        ok = True
        for v in exit_edge:
            seen, u = {v}, v
            while u != last:
                u = exit_edge[u]
                if u in seen:
                    ok = False
                    break
                seen.add(u)
            if not ok:
                break
        if ok:
            break
    for v, nxt in succ.items():
        if v != last:
            nxt.remove(exit_edge[v])
        random.shuffle(nxt)
        if v != last:
            nxt.append(exit_edge[v])
    out, v = [seq[0]], seq[0]
    for _ in range(len(seq) - 1):
        v = succ[v].pop(0)
        out.append(v)
    return ''.join(out)

def shuffle_seq(seq, preserve='mono'):
    '''Shuffle a sequence preserving the chosen statistic ('mono' or 'di').'''
    if preserve == 'mono':
        chars = list(seq)
        random.shuffle(chars)
        return ''.join(chars)
    if preserve == 'di':
        return dinuc_shuffle(seq)
    raise ValueError("preserve must be 'mono' or 'di'")

def empirical_pvalue(seq1, seq2, aligner, n_shuffles=1000, seed=42, preserve='mono'):
    '''Compute empirical p-value from a shuffled null distribution.

    Returns observed score, p-value, and the full null distribution.
    preserve='di' keeps dinucleotide counts (use for DNA).
    The (n_at_or_above + 1) / (n_shuffles + 1) correction (Phipson & Smyth 2010)
    avoids reporting p == 0 when no shuffle exceeds the observed score.
    '''
    random.seed(seed)
    observed = aligner.score(seq1, seq2)
    null_scores = [aligner.score(shuffle_seq(seq1, preserve), seq2) for _ in range(n_shuffles)]
    n_at_or_above = sum(1 for s in null_scores if s >= observed)
    pvalue = (n_at_or_above + 1) / (n_shuffles + 1)
    return observed, pvalue, null_scores

if __name__ == '__main__':
    aligner = PairwiseAligner(mode='local',
                              substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-11, extend_gap_score=-1)
    observed, p, _ = empirical_pvalue('MKTIIALSYIFCLVFA', 'MKAIIVCSCLLVFFA', aligner, n_shuffles=1000)
    print(f'Observed: {observed:.1f}  Empirical p-value: {p:.4f}')
