'''Match query MS/MS against a reference library with a score floor, a matched-peak floor and tie detection.

Inputs:  REFERENCES.mgf  library spectra (each needs a compound_name and a precursor m/z)
         QUERIES.mgf     query spectra
Usage:   python scripts/match_library.py REFERENCES.mgf QUERIES.mgf [--tolerance 0.005]
Prints one line per query: Level 2a candidate (single hit), "tied -> Level 3 (isomer wall)"
(several hits within TIE_MARGIN of the top score) or "no confident candidate -> Level 5".
Checked on matchms 0.33.1.
'''
import argparse

from matchms import calculate_scores
from matchms.filtering import default_filters, normalize_intensities, add_precursor_mz
from matchms.importing import load_from_mgf
try:
    from matchms.similarity import ModifiedCosineGreedy as ModifiedCosine  # matchms 0.33+
except ImportError:
    from matchms.similarity import ModifiedCosine          # matchms <= 0.32

SCORE_FLOOR = 0.7   # GNPS defaults (Wang 2016)
MATCH_FLOOR = 6
TIE_MARGIN = 0.02   # candidates within this margin of the top score are tied, not resolved


def prepare(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = add_precursor_mz(spectrum)  # required for ModifiedCosine; it only logs a warning
                                            # if precursor_mz is not derivable -- the
                                            # AssertionError fires later, in calculate_scores()
    return normalize_intensities(spectrum)


def match(references_raw, queries_raw, tolerance=0.005):
    queries = [prepare(s) for s in queries_raw]
    references = [prepare(s) for s in references_raw]
    scores = calculate_scores(references, queries, ModifiedCosine(tolerance=tolerance))
    calls = {}
    # CosineGreedy/ModifiedCosine return a structured array; the field names are
    # class-prefixed and version-dependent (e.g. 'ModifiedCosineGreedy_score' in 0.33),
    # so derive them from the dtype rather than hard-coding.
    for query in queries:
        pairs = scores.scores_by_query(query)
        name = query.get('compound_name')
        if not pairs:  # no reference shares a peak with this query: matchms stores no scores at all
            print(name, "-> no confident candidate -> Level 5")
            calls[name] = ('5', [])
            continue
        score_field, match_field = pairs[0][1].dtype.names
        ranked = sorted(pairs, key=lambda pair: pair[1][score_field], reverse=True)
        top_score = ranked[0][1][score_field]
        # Keep every hit within TIE_MARGIN of the top score, not just the argmax -- isomers
        # routinely score identically (the isomer wall), and taking a single winner
        # would silently launder a tie into a false Level 2a identification.
        passing = [(ref, hit) for ref, hit in ranked
                   if hit[score_field] >= top_score - TIE_MARGIN
                   and hit[score_field] >= SCORE_FLOOR and hit[match_field] >= MATCH_FLOOR]
        if len(passing) > 1:
            print([ref.get('compound_name') for ref, _ in passing], "tied -> Level 3 (isomer wall)")
            calls[name] = ('3', [ref.get('compound_name') for ref, _ in passing])
        elif len(passing) == 1:
            ref, hit = passing[0]
            print(ref.get('compound_name'), hit[score_field], hit[match_field])  # Level 2a candidate
            calls[name] = ('2a', [ref.get('compound_name')])
        else:
            print(name, "-> no confident candidate -> Level 5")
            calls[name] = ('5', [])
    return calls


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('references')
    ap.add_argument('queries')
    ap.add_argument('--tolerance', type=float, default=0.005)
    args = ap.parse_args()
    match(list(load_from_mgf(args.references)), list(load_from_mgf(args.queries)), args.tolerance)
