'''Bayesian protein inference with EPIFANY (pyOpenMS BayesianProteinInferenceAlgorithm).

Purpose: propagate belief over the peptide-protein graph and print each indistinguishable group with
its posterior. Input idXML must already carry PSM posterior error probabilities (Percolator or
IDPosteriorErrorProbability). The TOPP tool of the same method is `Epifany`.
Usage: python epifany_inference.py peptides_with_pep.idXML [--greedy-group-resolution] [--out groups.tsv]
Checked on pyOpenMS 3.5.0.
'''
import argparse

from pyopenms import IdXMLFile, BayesianProteinInferenceAlgorithm, PeptideIdentificationList

ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('idxml')
ap.add_argument('--greedy-group-resolution', action='store_true',
                help='third positional argument of inferPosteriorProbabilities; changed nothing on the reference idXML')
ap.add_argument('--out', help='optional TSV: accessions (;-joined) and group probability')
args = ap.parse_args()

protein_ids = []
peptide_ids = PeptideIdentificationList()
IdXMLFile().load(args.idxml, protein_ids, peptide_ids)

algo = BayesianProteinInferenceAlgorithm()
# EPIFANY expects PSM posteriors as input. The third POSITIONAL argument is
# greedy_group_resolution. Unlike BasicProteinInferenceAlgorithm's parameter of the same
# name, flipping it here did not change the result on the reference idXML (583 groups at 1%
# picked FDR, 3.43% true FDP with either value) -- so do NOT assume it removed subsumable
# proteins. Check the surviving groups, or use Basic + greedy when FDP is what you care about.
algo.inferPosteriorProbabilities(protein_ids, peptide_ids, args.greedy_group_resolution)

rows = []
for group in protein_ids[0].getIndistinguishableProteins():
    accs = [a.decode() for a in group.accessions]
    rows.append((accs, group.probability))
    print(accs, group.probability)   # higher = more likely present
if args.out:
    with open(args.out, 'w', encoding='utf-8') as fh:
        fh.write('accessions\tprobability\n')
        for accs, p in rows:
            fh.write(';'.join(accs) + f'\t{p}\n')
