'''Protein groups from an FDR-filtered peptide idXML: score aggregation + greedy resolution + picked-group FDR.

Purpose: run pyOpenMS BasicProteinInferenceAlgorithm with indistinguishable-group annotation AND
greedy group resolution, apply picked protein-group FDR, and emit one record per group:
leading_protein, accessions, n_peptides, n_unique_peptides (unique to the GROUP), is_decoy, qvalue.
Inputs: idXML whose PSMs carry PEPs and that KEEPS the decoy hits; the decoy prefix of the search
(DECOY_ OpenMS/Comet, rev_ Philosopher, REV__ MaxQuant).
Usage: python group_proteins_pyopenms.py peptides_1pct_fdr.idXML [--decoy-prefix DECOY_] [--fdr 0.01] [--out groups.tsv]
Checked on pyOpenMS 3.5.0.
'''
import argparse
import csv

from pyopenms import (IdXMLFile, BasicProteinInferenceAlgorithm, PeptideIdentificationList,
                      FalseDiscoveryRate, String)

ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('idxml')
ap.add_argument('--decoy-prefix', default='DECOY_', help='tool-specific: DECOY_ OpenMS/Comet, rev_ Philosopher, REV__ MaxQuant')
ap.add_argument('--fdr', type=float, default=0.01)
ap.add_argument('--out', help='optional TSV of every group (targets and decoys)')
args = ap.parse_args()
DECOY_PREFIX = args.decoy_prefix

protein_ids = []
peptide_ids = PeptideIdentificationList()   # pyOpenMS 3.5+: a plain [] fails
# protein_ids is FIRST in both load() and store() for IdXMLFile
IdXMLFile().load(args.idxml, protein_ids, peptide_ids)

inference = BasicProteinInferenceAlgorithm()
params = inference.getParameters()
# annotate_indistinguishable_groups reports indistinguishable proteins as ONE group
params.setValue('annotate_indistinguishable_groups', 'true')
# greedy_group_resolution assigns shared peptides to the best group, so subsumable proteins drop out
params.setValue('greedy_group_resolution', 'true')
inference.setParameters(params)
inference.run(peptide_ids, protein_ids)

# picked protein-group FDR: (decoy string, is prefix, groups too); group.probability becomes a q-value
FalseDiscoveryRate().applyPickedProteinFDR(protein_ids[0], String(DECOY_PREFIX), True, True)

# peptide sequences behind each protein, from the best hit of every PSM
peps_of = {}
for pid in peptide_ids:
    hit = pid.getHits()[0]
    for ev in hit.getPeptideEvidences():
        peps_of.setdefault(ev.getProteinAccession(), set()).add(hit.getSequence().toString())

# group record: leading_protein, accessions, n_peptides, n_unique_peptides, is_decoy, qvalue
# (the same record examples/protein_groups.py returns)
groups = []
for g in protein_ids[0].getIndistinguishableProteins():
    accs = [a.decode() for a in g.accessions]   # bytes; pyOpenMS sorts them alphabetically
    peptides = set().union(*(peps_of.get(a, set()) for a in accs))
    others = set().union(*(s for a, s in peps_of.items() if a not in accs))   # unique = not in any other protein
    groups.append({
        # members share the same evidence; choose the lead explicitly: canonical (no -N isoform suffix) first
        'leading_protein': sorted(accs, key=lambda a: ('-' in a, a))[0],
        'accessions': accs,
        'n_peptides': len(peptides),
        'n_unique_peptides': len(peptides - others),
        'is_decoy': all(a.startswith(DECOY_PREFIX) for a in accs),
        'qvalue': g.probability,
    })
passing = [g for g in groups if not g['is_decoy'] and g['qvalue'] <= args.fdr]
print(len(passing), f'groups at {args.fdr:.0%} picked-group FDR')
for g in passing[:5]:
    print(g)

if args.out:
    with open(args.out, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['leading_protein', 'accessions', 'n_peptides', 'n_unique_peptides', 'is_decoy', 'qvalue'])
        for g in groups:
            w.writerow([g['leading_protein'], ';'.join(g['accessions']), g['n_peptides'],
                        g['n_unique_peptides'], g['is_decoy'], g['qvalue']])
