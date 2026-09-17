'''Infinite-sites plot: regress each dated node's posterior HPD width against its posterior
mean age, parsed straight from one or more MCMCTree out.txt files (dos Reis and Yang 2011).

Under infinite sequence data the relationship becomes linear THROUGH THE ORIGIN, with the
residual width set entirely by calibration uncertainty -- see SKILL.md's "More sequence data
does not escape the calibration" rule. This script fits that line (forced through the origin,
as the theory predicts) to one out.txt and reports the slope and R^2; passed two or more
out.txt files from runs with increasing data (e.g. 1 locus, then 2 loci, same tree and
calibrations), it also reports whether the slope has stabilized and R^2 is rising toward 1 --
that combination means the points are settling onto the calibration-limited line, so more
sequence data will not narrow these dates further. A single run with high R^2 is suggestive but
not conclusive on its own; the trend across runs is the real test.

Usage:
  python infinite_sites_plot.py post/out.txt [post2/out.txt ...]
    out.txt   MCMCTree's output file (usedata = 2 in.BV posterior run); needs >= 3 dated
              (t_n*) nodes with a "Posterior means (95% Equal-tail CI) (95% HPD CI)" table.
'''
# Reference: PAML/MCMCTree 4.10.10 | Verify out.txt's column layout if version differs

import re
import sys

# t_nN  mean  ( eq-tail lo,  eq-tail hi) ( hpd lo,  hpd hi)  [hpd width]
NODE_ROW = re.compile(
    r'^(t_n\d+)\s+([\d.eE+-]+)\s+\(\s*([\d.eE+-]+),\s*([\d.eE+-]+)\)\s+\(\s*([\d.eE+-]+),\s*([\d.eE+-]+)\)',
    re.M)


def parse_out_txt(path):
    '''Return {node: (posterior_mean, hpd_lo, hpd_hi)} from MCMCTree's out.txt summary table.'''
    text = open(path, encoding='utf-8', errors='replace').read()
    nodes = {m.group(1): (float(m.group(2)), float(m.group(5)), float(m.group(6)))
              for m in NODE_ROW.finditer(text)}
    if len(nodes) < 3:
        raise ValueError(f'{path}: only {len(nodes)} dated (t_n*) nodes found; need >= 3')
    return nodes


def infinite_sites_fit(nodes):
    '''Slope and R^2 of HPD width on posterior mean age, forced through the origin
    (dos Reis and Yang 2011): slope = sum(mean*width) / sum(mean**2); R^2 measured against
    the mean-width baseline, same as a normal regression's R^2.'''
    means = [v[0] for v in nodes.values()]
    widths = [hi - lo for _, lo, hi in nodes.values()]
    sum_mm = sum(x * x for x in means)
    if sum_mm == 0:
        raise ValueError('all posterior means are zero')
    slope = sum(x * y for x, y in zip(means, widths)) / sum_mm
    mean_w = sum(widths) / len(widths)
    ss_tot = sum((w - mean_w) ** 2 for w in widths)
    ss_res = sum((w - slope * x) ** 2 for x, w in zip(means, widths))
    r2 = 1 - ss_res / ss_tot if ss_tot else float('nan')
    return slope, r2, len(means)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    print(f"{'out.txt':45s} {'n_nodes':>7s} {'slope (width/age)':>18s} {'R^2':>8s}")
    fits = []
    for path in sys.argv[1:]:
        slope, r2, n = infinite_sites_fit(parse_out_txt(path))
        fits.append((path, slope, r2))
        print(f'{path:45s} {n:7d} {slope:18.4f} {r2:8.4f}')
    if len(fits) >= 2:
        _, slope0, r2_0 = fits[0]
        _, slope1, r2_1 = fits[-1]
        slope_stable = abs(slope1 - slope0) < 0.25 * max(abs(slope0), 1e-9)
        r2_rising = r2_1 >= r2_0
        if slope_stable and r2_rising:
            print('Slope stable and R^2 rising toward 1 across runs: points are settling onto the')
            print('through-origin line. More sequence data will NOT narrow these dates further --')
            print('the residual HPD width is set by calibration uncertainty. Invest in better')
            print('fossils, not more sites.')
        else:
            print('Slope/R^2 have not stabilized across these runs -- more loci may still narrow')
            print('the dates; re-run this diagnostic as data accumulate.')
    else:
        _, slope, r2 = fits[0]
        print(f'Single run: slope {slope:.4f}, R^2 {r2:.4f}. R^2 alone from one run is suggestive,')
        print('not conclusive -- compare against a run with more loci to see whether the slope has')
        print('stabilized before concluding you are at the calibration floor.')
