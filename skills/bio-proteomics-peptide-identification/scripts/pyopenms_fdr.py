"""Annotate target/decoy and estimate FDR with pyOpenMS, keeping PSMs at q <= 0.01.

Input : idXML from pyopenms_search.py (or another engine), and the same target+decoy FASTA.
Output: idXML of target PSMs at the FDR threshold (q-values as the score).
Usage : python pyopenms_fdr.py search_results.idXML human_target_decoy.fasta psms_1pct.idXML
        [--decoy-string DECOY_] [--fdr 0.01]
Checked on pyOpenMS 3.5.0.
"""
import argparse

from pyopenms import PeptideIndexing, FalseDiscoveryRate, IDFilter, FASTAFile, IdXMLFile, PeptideIdentificationList

ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('in_idxml')
ap.add_argument('fasta', help='the target + decoy FASTA the search used')
ap.add_argument('out_idxml')
ap.add_argument('--decoy-string', default='DECOY_', help='must match the decoy prefix in the FASTA')
ap.add_argument('--fdr', type=float, default=0.01)
args = ap.parse_args()

protein_ids = []
peptide_ids = PeptideIdentificationList()
IdXMLFile().load(args.in_idxml, protein_ids, peptide_ids)   # protein_ids FIRST
n_in = peptide_ids.size()

fasta = []
FASTAFile().load(args.fasta, fasta)
indexer = PeptideIndexing()
params = indexer.getParameters()
params.setValue('decoy_string', args.decoy_string)
params.setValue('decoy_string_position', 'prefix')
indexer.setParameters(params)
indexer.run(fasta, protein_ids, peptide_ids)   # sets target/decoy flags on every hit

FalseDiscoveryRate().apply(peptide_ids)         # concatenated competition -> per-PSM q-value as the new score
IDFilter().filterHitsByScore(peptide_ids, args.fdr)   # 0.01 = 1% FDR, the community list-level standard
IDFilter().removeDecoyHits(peptide_ids)
IDFilter().removeEmptyIdentifications(peptide_ids)   # spectra whose only hits were filtered out

IdXMLFile().store(args.out_idxml, protein_ids, peptide_ids)
print(f'{n_in} PSMs in, {peptide_ids.size()} target PSMs at q <= {args.fdr} -> {args.out_idxml}')
