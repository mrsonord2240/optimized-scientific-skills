'''Date-randomization test (Duchene et al. 2015), the formal temporal-signal test.

A positive root-to-tip slope (root_to_tip.py) is necessary but not sufficient:
short-time-span data can still give a positive slope from root placement alone.
This shuffles tip sampling dates N times, reruns LSD2 (via IQ-TREE2 --date) on
the SAME fixed tree topology each time, and checks whether the real-data rate
estimate falls outside the range of shuffled-replicate rates. If the real rate
sits inside that cloud, there is no temporal signal and any date from that
dataset is an artifact of the clock-rate prior, not the data.

Usage:
  python date_randomization.py ALN.FA TREEFILE DATES.TSV [N_REPLICATES] [SEED]
    ALN.FA      alignment IQ-TREE used to build TREEFILE
    TREEFILE    Newick tree, topology held fixed across all replicates (-te)
    DATES.TSV   tip <TAB> decimal-year-date (same file used for the real --date run)
    N_REPLICATES  default 20
    SEED          default 42

Requires `iqtree2` (or `iqtree3` on IQ-TREE 3.x) on PATH, or edit IQTREE below.
'''
# Reference: IQ-TREE 2.4.0 / LSD2 2.4.4 | Verify --date flags if version differs

import random
import re
import shutil
import subprocess
import sys

IQTREE = shutil.which('iqtree2') or shutil.which('iqtree3') or shutil.which('iqtree') or 'iqtree2'
RATE_RE = re.compile(r'rate ([\d.eE+-]+), tMRCA ([\d.eE+-]+)')


def load_dates(path):
    tips, vals = [], []
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            tip, date = line.split('\t')
            tips.append(tip)
            vals.append(date)
    return tips, vals


def write_dates(path, tips, vals):
    with open(path, 'w', encoding='utf-8', newline='\n') as fh:
        for t, v in zip(tips, vals):
            fh.write(f'{t}\t{v}\n')


def run_lsd(aln, tree, dates_path, prefix):
    '''Run LSD2 (via IQ-TREE2 --date) on a FIXED tree (-te) and parse the rate estimate.'''
    cmd = [IQTREE, '-s', aln, '-m', 'GTR+G', '-te', tree, '--date', dates_path,
           '--prefix', prefix, '-redo', '-T', '1']
    r = subprocess.run(cmd, capture_output=True, text=True)
    m = RATE_RE.search(r.stdout)
    if not m:
        raise RuntimeError('no rate parsed from iqtree2 output:\n' + r.stdout[-1500:] + r.stderr[-1500:])
    return float(m.group(1)), float(m.group(2))


def randomization_test(aln, tree, dates_path, n_rep, seed):
    rng = random.Random(seed)
    tips, vals = load_dates(dates_path)
    real_rate, real_tmrca = run_lsd(aln, tree, dates_path, 'randtest_real')
    shuf_rates = []
    for i in range(n_rep):
        shuffled = vals[:]
        rng.shuffle(shuffled)
        dpath = f'randtest_shuf{i}.tsv'
        write_dates(dpath, tips, shuffled)
        rate, _ = run_lsd(aln, tree, dpath, f'randtest_shuf{i}')
        shuf_rates.append(rate)
    return real_rate, real_tmrca, shuf_rates


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    aln, tree, dates_path = sys.argv[1], sys.argv[2], sys.argv[3]
    n_rep = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 42

    real_rate, real_tmrca, shuf_rates = randomization_test(aln, tree, dates_path, n_rep, seed)
    print(f'real: rate={real_rate:.6e} tMRCA={real_tmrca:.3f}')
    for i, r in enumerate(shuf_rates):
        print(f'shuffle {i}: rate={r:.6e}')

    lo, hi = min(shuf_rates), max(shuf_rates)
    outside = real_rate < lo or real_rate > hi
    print(f'\nshuffled rate range over {n_rep} replicates: [{lo:.6e}, {hi:.6e}]')
    print(f"real rate {real_rate:.6e} is {'OUTSIDE' if outside else 'INSIDE'} the randomized range "
          f"-> temporal signal {'CONFIRMED -- safe to date' if outside else 'NOT CONFIRMED -- do not report a date'}")
