"""Inspect mzML/mzXML spectra without treating missing precursor metadata as zero.

Usage: python inspect_mzml.py sample.mzML-or-mzXML
"""
import sys
from pathlib import Path
from pyopenms import MSExperiment, MzMLFile, MzXMLFile


def inspect_mzml(path):
    exp = MSExperiment()
    suffix = Path(path).suffix.lower()
    if suffix == '.mzml':
        loader = MzMLFile()
    elif suffix == '.mzxml':
        loader = MzXMLFile()
    else:
        raise ValueError('Expected an .mzML or .mzXML file')
    loader.load(str(path), exp)
    counts = {'ms1': 0, 'ms2': 0, 'ms2_without_precursor': 0, 'offsets_unknown': 0}
    for spectrum in exp:
        if spectrum.getMSLevel() == 1:
            counts['ms1'] += 1
            continue
        if spectrum.getMSLevel() != 2:
            continue
        counts['ms2'] += 1
        precursors = spectrum.getPrecursors()
        if not precursors:
            counts['ms2_without_precursor'] += 1
            continue
        precursor = precursors[0]
        width = precursor.getIsolationWindowLowerOffset() + precursor.getIsolationWindowUpperOffset()
        if width == 0:
            counts['offsets_unknown'] += 1
    return counts


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python inspect_mzml.py sample.mzML-or-mzXML')
    print(inspect_mzml(sys.argv[1]))
