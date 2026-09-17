# Reference: biopython 1.83+, pandas 2.2+, numpy 1.26+ | Verify API if version differs
import re
import numpy as np
import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

# === CONFIGURATION ===
output_dir = Path('library_design/')
output_dir.mkdir(exist_ok=True)

np.random.seed(42)

# === SKILL.md's documented design functions (verbatim; see SKILL.md's
# "Score and Rank sgRNAs for a Target Gene" section for the source of truth) ===

def find_sgrna_candidates(cds_sequence, pam='NGG', guide_length=20):
    '''Return all protospacer candidates with PAM coordinates on + strand.
    Caller must filter by exon position and Azimuth/CFD score.'''
    pam_pattern = re.compile(f'(?=([ACGT]{{{guide_length}}}{pam.replace("N", "[ACGT]")}))')
    candidates = []
    for strand, seq in [('+', cds_sequence), ('-', str(Seq(cds_sequence).reverse_complement()))]:
        for m in pam_pattern.finditer(seq):
            spacer = m.group(1)[:guide_length]
            if 'TTTT' in spacer or spacer.count('G') + spacer.count('C') not in range(6, 15):
                continue
            candidates.append({'spacer': spacer, 'strand': strand,
                               'pos_in_cds': m.start() if strand == '+' else len(seq) - m.start() - 23,
                               'gc_frac': (spacer.count('G') + spacer.count('C')) / guide_length})
    return pd.DataFrame(candidates)

def annotate_exon_position(candidates_df, cds_length):
    '''Filter to protospacers within first 5-65% of CDS (Brunello convention).'''
    lo, hi = 0.05 * cds_length, 0.65 * cds_length
    return candidates_df[(candidates_df['pos_in_cds'] >= lo) & (candidates_df['pos_in_cds'] <= hi)].copy()

def select_independent_guides(candidates_df, n_guides, min_spacing=5, score_col='score'):
    '''Greedily pick up to n_guides candidates that are >=min_spacing nt apart
    (see SKILL.md for why composition filters alone let near-duplicates through).'''
    ranked = candidates_df.sort_values(score_col, ascending=False)
    selected = []
    for _, cand in ranked.iterrows():
        if any(abs(cand['pos_in_cds'] - s['pos_in_cds']) < min_spacing for s in selected):
            continue
        selected.append(cand)
        if len(selected) == n_guides:
            break
    return pd.DataFrame(selected)

def generate_synthetic_cds(length=600):
    '''Demo-only stand-in for a real Ensembl/RefSeq CDS pull. A real design
    run replaces this with the gene's actual coding sequence.'''
    bases = ['A', 'C', 'G', 'T']
    return ''.join(np.random.choice(bases, length))

def generate_sgrna_sequence(length=20):
    bases = ['A', 'C', 'G', 'T']
    return ''.join(np.random.choice(bases, length))

def score_sgrna(sequence):
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence)
    gc_score = 1 - abs(gc_content - 0.5) * 2
    poly_t_penalty = 0 if 'TTTT' in sequence else 1
    start_g_bonus = 1 if sequence.startswith('G') else 0.8
    return gc_score * poly_t_penalty * start_g_bonus, gc_content

# === 1. DEFINE TARGET GENES ===
print('Defining target gene list...')

target_genes = ['TP53', 'BRCA1', 'BRCA2', 'KRAS', 'NRAS', 'BRAF', 'MYC', 'MYCN', 'CDK4', 'CDK6',
                'RB1', 'PTEN', 'PIK3CA', 'AKT1', 'MTOR', 'EGFR', 'ERBB2', 'MET', 'ALK', 'ROS1']

print(f'Target genes: {len(target_genes)}')

# === 2. DESIGN sgRNAs USING THE DOCUMENTED CDS/EXON-POSITION METHOD ===
# (previously this demo generated fully random 20nt sequences and scored them
# by GC content alone, which never exercised find_sgrna_candidates /
# annotate_exon_position / select_independent_guides -- the actual method
# documented in SKILL.md. It now calls that method against a synthetic CDS
# per gene, same as SKILL.md's own worked examples.)
print('\nDesigning sgRNAs for each gene from a synthetic CDS...')

guides_per_gene = 4
all_guides = []
shortfall_genes = []

for gene in target_genes:
    cds = generate_synthetic_cds(length=600)
    candidates = find_sgrna_candidates(cds)
    candidates = annotate_exon_position(candidates, len(cds))
    if candidates.empty:
        shortfall_genes.append((gene, 0))
        continue
    candidates['score'] = 1 - (candidates['gc_frac'] - 0.5).abs() * 2
    selected = select_independent_guides(candidates, guides_per_gene, min_spacing=5, score_col='score')
    if len(selected) < guides_per_gene:
        shortfall_genes.append((gene, len(selected)))

    for i, (_, guide) in enumerate(selected.iterrows()):
        all_guides.append({
            'gene': gene,
            'guide_number': i + 1,
            'sequence': guide['spacer'],
            'score': guide['score'],
            'gc_content': guide['gc_frac'],
            'type': 'targeting'
        })

print(f'Targeting guides: {len(all_guides)}')
if shortfall_genes:
    print(f'Genes below the {guides_per_gene}-guide quota after independence filtering: {shortfall_genes}')

# === 3. ADD CONTROLS ===
# NOTE: absolute control counts below are sized for this 20-gene demo only.
# At genome scale, scale non-targeting controls to ~1% of the full library
# (500-1,000 in a 70k-guide library) and safe-harbor/essential/non-essential
# controls to 50-100 each -- see SKILL.md's "Control Guides" table. Do not
# read this demo's ~44% control fraction as a real-library target.
print('\nAdding control guides...')

n_nontargeting = 50
for i in range(n_nontargeting):
    while True:
        seq = generate_sgrna_sequence()
        score, gc = score_sgrna(seq)
        if score > 0.5:
            break

    all_guides.append({
        'gene': f'NonTargeting_{i+1:03d}',
        'guide_number': 1,
        'sequence': seq,
        'score': 0,
        'gc_content': gc,
        'type': 'non-targeting'
    })

essential_genes = ['RPL11', 'RPS3', 'EIF3A', 'POLR2A', 'CDK1', 'SF3B1', 'U2AF1', 'PRPF8', 'SNRPD1', 'SNRPE']
for gene in essential_genes:
    seq = generate_sgrna_sequence()
    while 'TTTT' in seq:                     # poly-T terminates U6; same filter as targeting guides
        seq = generate_sgrna_sequence()
    _, gc = score_sgrna(seq)
    all_guides.append({
        'gene': gene,
        'guide_number': 1,
        'sequence': seq,
        'score': 1.0,
        'gc_content': gc,
        'type': 'essential-control'
    })

for gene in ['AAVS1', 'ROSA26']:
    seq = generate_sgrna_sequence()
    while 'TTTT' in seq:                     # poly-T terminates U6; same filter as targeting guides
        seq = generate_sgrna_sequence()
    _, gc = score_sgrna(seq)
    all_guides.append({
        'gene': gene,
        'guide_number': 1,
        'sequence': seq,
        'score': 1.0,
        'gc_content': gc,
        'type': 'safe-harbor'
    })

library_df = pd.DataFrame(all_guides)
print(f'Total library size: {len(library_df)}')

# === 4. DESIGN OLIGOS ===
print('\nDesigning cloning oligos...')

forward_prefix = 'CACCG'
reverse_prefix = 'AAAC'
reverse_suffix = 'C'

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

library_df['forward_oligo'] = forward_prefix + library_df['sequence']
library_df['reverse_oligo'] = reverse_prefix + library_df['sequence'].apply(reverse_complement) + reverse_suffix

# === 5. LIBRARY QC ===
print('\n=== LIBRARY QC ===')

print(f"\nGuide type distribution:")
print(library_df['type'].value_counts())

targeting = library_df[library_df['type'] == 'targeting']
print(f"\nTargeting guides per gene:")
guides_per = targeting.groupby('gene').size()
print(f"  Mean: {guides_per.mean():.1f}")
print(f"  Min: {guides_per.min()}")
print(f"  Max: {guides_per.max()}")

print(f"\nGC content distribution:")
print(f"  Mean: {library_df['gc_content'].mean():.1%}")
print(f"  Std: {library_df['gc_content'].std():.1%}")
print(f"  Range: {library_df['gc_content'].min():.1%} - {library_df['gc_content'].max():.1%}")

poly_t_count = library_df['sequence'].apply(lambda x: 'TTTT' in x).sum()
print(f"\nPoly-T sequences: {poly_t_count} ({poly_t_count/len(library_df):.1%})")

control_pct = (library_df['type'] != 'targeting').sum() / len(library_df)
print(f"Control fraction: {control_pct:.1%} (demo-scale only -- see NOTE above; real libraries target ~1% NTC)")

# === 6. EXPORT ===
print('\n=== EXPORTING LIBRARY ===')

library_df.to_csv(output_dir / 'library_design.csv', index=False)

oligo_order = library_df[['gene', 'guide_number', 'forward_oligo', 'reverse_oligo', 'type']].copy()
oligo_order['oligo_id'] = library_df['gene'] + '_g' + library_df['guide_number'].astype(str)
oligo_order.to_csv(output_dir / 'oligo_order.csv', index=False)

summary = {
    'total_guides': len(library_df),
    'targeting_genes': library_df[library_df['type'] == 'targeting']['gene'].nunique(),
    'targeting_guides': (library_df['type'] == 'targeting').sum(),
    'nontargeting_controls': (library_df['type'] == 'non-targeting').sum(),
    'essential_controls': (library_df['type'] == 'essential-control').sum(),
    'safeharbor_controls': (library_df['type'] == 'safe-harbor').sum(),
    'mean_gc': library_df['gc_content'].mean(),
    'guides_per_gene': guides_per_gene
}

pd.DataFrame([summary]).to_csv(output_dir / 'library_summary.csv', index=False)

print(f"\nLibrary design saved to {output_dir}/")
print(f"  - library_design.csv: Full library with scores")
print(f"  - oligo_order.csv: Oligo ordering format")
print(f"  - library_summary.csv: Summary statistics")
