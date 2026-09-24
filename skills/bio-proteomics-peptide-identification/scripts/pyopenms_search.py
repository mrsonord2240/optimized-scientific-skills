"""Database search with pyOpenMS SimpleSearchEngineAlgorithm -> scored PSMs as idXML.

Input : centroided mzML, and a FASTA that already holds target + decoy sequences
        (decoys carry a recognisable prefix, e.g. DECOY_).
Output: idXML with the target_decoy flag on every hit (feed it to pyopenms_fdr.py).
Usage : python pyopenms_search.py sample.mzML human_target_decoy.fasta search_results.idXML
        [--precursor-ppm 10] [--fragment-da 0.02] [--missed-cleavages 2]
Checked on pyOpenMS 3.5.0.
"""
import argparse

from pyopenms import SimpleSearchEngineAlgorithm, IdXMLFile, PeptideIdentificationList

ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
ap.add_argument('mzml')
ap.add_argument('fasta', help='target + decoy FASTA')
ap.add_argument('out_idxml')
ap.add_argument('--precursor-ppm', type=float, default=10.0)
ap.add_argument('--fragment-da', type=float, default=0.02, help='HCD Orbitrap; the engine default is 10 ppm')
ap.add_argument('--missed-cleavages', type=int, default=2, help='the engine default is 1')
args = ap.parse_args()

protein_ids = []
peptide_ids = PeptideIdentificationList()   # pyOpenMS 3.5+: a plain [] raises TypeError
search = SimpleSearchEngineAlgorithm()
p = search.getParameters()
p.setValue('precursor:mass_tolerance', args.precursor_ppm)
p.setValue('precursor:mass_tolerance_unit', 'ppm')
p.setValue('fragment:mass_tolerance', args.fragment_da)
p.setValue('fragment:mass_tolerance_unit', 'Da')
p.setValue('peptide:missed_cleavages', args.missed_cleavages)
p.setValue('modifications:fixed', [b'Carbamidomethyl (C)'])
p.setValue('modifications:variable', [b'Oxidation (M)'])
search.setParameters(p)
# spectra are scored against in-silico fragment ions of every candidate peptide
search.search(args.mzml, args.fasta, protein_ids, peptide_ids)

# protein_ids FIRST in load/store -- the OpenMS argument order is fixed
IdXMLFile().store(args.out_idxml, protein_ids, peptide_ids)
print(f'{peptide_ids.size()} PSMs written to {args.out_idxml}')
