'''PAML codeml ancestral reconstruction with site-confidence parsing'''

import subprocess
import re
import os
from collections import defaultdict
from paml_rst import parse_marginal_rst


def write_codeml_ctl(alignment, tree, out_dir, seqtype='codon', model='M0'):
    '''Write codeml control file for ancestral reconstruction.

    seqtype: 1=codon, 2=aa, 3=codon translated
    For protein resurrection use seqtype=2 with model=3 (empirical+gamma).
    RateAncestor=1 produces rst file with per-site posteriors.
    '''
    if seqtype == 'codon':
        ctl = f'''
      seqfile = {alignment}
     treefile = {tree}
      outfile = {out_dir}/mlc
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
 RateAncestor = 1
    cleandata = 0
   getSE = 1
'''
    else:
        ctl = f'''
      seqfile = {alignment}
     treefile = {tree}
      outfile = {out_dir}/mlc
      runmode = 0
      seqtype = 2
        model = 3
   aaRatefile = lg.dat
       alpha = 0.5
        ncatG = 4
 RateAncestor = 1
    cleandata = 0
'''
    ctl_path = os.path.join(out_dir, 'codeml.ctl')
    open(ctl_path, 'w').write(ctl)
    return ctl_path


def parse_rst_posteriors(rst_file):
    """Return marginal best-state probabilities grouped by ancestral node."""
    _, parsed = parse_marginal_rst(rst_file)
    result = {}
    for node, rows in parsed.items():
        result[int(node.removeprefix('Node'))] = []
        for row in rows:
            item = {'site': row['site'], 'state': row['state'], 'prob': row['probability']}
            if 'amino_acid_state' in row:
                item['codon_translation'] = row['codon_translation']
                item['amino_acid_state'] = row['amino_acid_state']
                item['amino_acid_probability'] = row['amino_acid_probability']
            result[int(node.removeprefix('Node'))].append(item)
    return result


def summarize_node_confidence(node_posteriors, ambig_cutoff=0.8):
    '''Per-node confidence summary for protein-resurrection construct design.'''
    summary = {}
    for node, sites in node_posteriors.items():
        probs = [s['prob'] for s in sites]
        n_high = sum(p >= 0.95 for p in probs)
        n_amb = sum(p < ambig_cutoff for p in probs)
        summary[node] = {
            'mean_post': sum(probs) / len(probs),
            'frac_high': n_high / len(probs),
            'n_ambiguous': n_amb,
            'ambiguous_sites': [s['site'] for s in sites if s['prob'] < ambig_cutoff]
        }
    return summary
