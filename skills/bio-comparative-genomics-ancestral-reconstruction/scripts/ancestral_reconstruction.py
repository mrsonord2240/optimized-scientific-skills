'''Ancestral sequence reconstruction'''
# Reference: biopython 1.83+, iq-tree 2.2+, paml 4.10+ | Verify API if version differs

import subprocess
import os
from paml_rst import parse_marginal_rst


def write_asr_control(alignment, tree, outfile, output_dir='.'):
    '''Generate PAML control file for ancestral reconstruction

    RateAncestor = 1 enables reconstruction
    Output goes to RST file in working directory
    '''
    ctl = f'''      seqfile = {alignment}
     treefile = {tree}
      outfile = {outfile}

        noisy = 3
      verbose = 1
      runmode = 0

      seqtype = 2
        model = 2
    aaRatefile = lg.dat

        clock = 0
        Mgene = 0

    fix_alpha = 0
        alpha = 0.5
       Malpha = 0
        ncatG = 4

 RateAncestor = 1
    cleandata = 0
'''
    ctl_path = os.path.join(output_dir, 'asr.ctl')
    with open(ctl_path, 'w') as f:
        f.write(ctl)

    return ctl_path


def parse_rst_ancestors(rst_file):
    """Parse reconstructed marginal MAP sequences from a current PAML rst."""
    ancestors, _ = parse_marginal_rst(rst_file)
    return ancestors


def extract_site_probabilities(rst_file):
    '''Extract posterior probabilities for ancestral states

    Confidence levels:
    - P > 0.95: High confidence (reliable for resurrection)
    - P 0.80-0.95: Moderate confidence
    - P < 0.80: Low confidence (consider alternatives)
    '''
    _, by_node = parse_marginal_rst(rst_file)
    return [record for records in by_node.values() for record in records]


def assess_reconstruction_quality(site_probs):
    '''Assess overall quality of ancestral reconstruction

    Returns quality metrics for experimental planning
    '''
    if not site_probs:
        return {
            'status': 'no_data', 'total_sites': 0, 'total_nodes': 0,
            'total_node_site_estimates': 0, 'high_confidence': 0,
            'moderate_confidence': 0, 'low_confidence': 0,
            'mean_probability': None, 'high_conf_fraction': None,
            'quality': 'unknown', 'message': 'No probability data found; reconstruction was not assessed.'
        }

    probs = [s['probability'] for s in site_probs]
    n_estimates = len(probs)
    total_sites = len({s.get('site') for s in site_probs if s.get('site') is not None})
    node_ids = {s.get('node') for s in site_probs if s.get('node') is not None}

    # Count sites by confidence
    high_conf = sum(1 for p in probs if p > 0.95)
    mod_conf = sum(1 for p in probs if 0.8 <= p <= 0.95)
    low_conf = sum(1 for p in probs if p < 0.8)

    mean_prob = sum(probs) / n_estimates
    high_conf_frac = high_conf / n_estimates

    # Overall quality assessment
    # For protein resurrection:
    # - >90% high confidence sites: Excellent, proceed with ML sequence
    # - 70-90% high confidence: Good, but test alternatives at uncertain sites
    # - <70% high confidence: Poor, may need better alignment/tree
    if high_conf_frac > 0.9:
        quality = 'excellent'
    elif high_conf_frac > 0.7:
        quality = 'good'
    elif high_conf_frac > 0.5:
        quality = 'moderate'
    else:
        quality = 'poor'

    return {
        'status': 'complete',
        'total_sites': total_sites or n_estimates,
        'total_nodes': len(node_ids) if node_ids else None,
        'total_node_site_estimates': n_estimates,
        'high_confidence': high_conf,
        'moderate_confidence': mod_conf,
        'low_confidence': low_conf,
        'mean_probability': mean_prob,
        'high_conf_fraction': high_conf_frac,
        'quality': quality
    }


def identify_ambiguous_sites(site_probs, threshold=0.8):
    '''Find sites with ambiguous ancestral states

    These sites should be tested with alternative states
    in resurrection experiments
    '''
    return [s for s in site_probs if s['probability'] < threshold]


def summarize_asr_results(ancestors, quality_metrics, ambiguous_sites):
    '''Print summary of ancestral reconstruction'''
    print('Ancestral Sequence Reconstruction Results')
    print('=' * 50)

    print(f"\nReconstructed {len(ancestors)} ancestral nodes")
    for node, seq in list(ancestors.items())[:3]:
        print(f"  {node}: {seq[:50]}..." if len(seq) > 50 else f"  {node}: {seq}")

    print(f"\nReconstruction status: {quality_metrics.get('status', 'unknown')}")
    print(f"Reconstruction Quality: {quality_metrics['quality'].upper()}")
    if quality_metrics.get('status') != 'complete':
        print(f"  {quality_metrics.get('message', 'No complete probability summary is available.')}")
        return
    print(f"  Nodes: {quality_metrics.get('total_nodes')}; site-node estimates: "
          f"{quality_metrics.get('total_node_site_estimates', quality_metrics['total_sites'])}")
    print(f"  Total sites: {quality_metrics['total_sites']}")
    print(f"  High confidence (P>0.95): {quality_metrics['high_confidence']} "
          f"({quality_metrics['high_conf_fraction']*100:.1f}%)")
    print(f"  Low confidence (P<0.80): {quality_metrics['low_confidence']}")
    print(f"  Mean probability: {quality_metrics['mean_probability']:.3f}")

    if ambiguous_sites:
        print(f"\nAmbiguous sites requiring attention: {len(ambiguous_sites)}")
        for site in ambiguous_sites[:5]:
            print(f"  Site {site['site']}: {site['state']} (P={site['probability']:.3f})")
        if len(ambiguous_sites) > 5:
            print(f"  ... and {len(ambiguous_sites) - 5} more")


if __name__ == '__main__':
    print('Ancestral Sequence Reconstruction')
    print('=' * 50)

    # Example results (simulated)
    example_ancestors = {
        'Node6': 'MKFLILLFNILCLFPVLAADYKDDDDKGENLYFQG',
        'Node7': 'MKFLILLFNILCLFPVLAADYKDDDDKGENLYFQG',
        'Node8': 'MKFLVLLFNILCLFPVLAADYKDDDDKGDNLYFQG',
    }

    example_probs = [
        {'site': 1, 'state': 'M', 'probability': 0.999},
        {'site': 2, 'state': 'K', 'probability': 0.987},
        {'site': 3, 'state': 'F', 'probability': 0.923},
        {'site': 4, 'state': 'L', 'probability': 0.756},  # Ambiguous
        {'site': 5, 'state': 'I', 'probability': 0.812},
        {'site': 6, 'state': 'L', 'probability': 0.634},  # Ambiguous
    ]

    quality = assess_reconstruction_quality(example_probs)
    ambiguous = identify_ambiguous_sites(example_probs, threshold=0.8)

    summarize_asr_results(example_ancestors, quality, ambiguous)

    print('\n\nTo run on real data:')
    print('1. Prepare protein MSA (PHYLIP format)')
    print('2. Generate phylogenetic tree')
    print('3. ctl = write_asr_control(aln, tree, "asr.mlc")')
    print('4. subprocess.run(["codeml", ctl])')
    print('5. ancestors = parse_rst_ancestors("rst")')
    print('6. probs = extract_site_probabilities("rst")')
