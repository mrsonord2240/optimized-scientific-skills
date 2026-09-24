'''Confirm ELink source UIDs resolve in the stated dbfrom, and show what they are.

ELink never errors on a wrong-namespace UID: it returns LinkSetDb == [] (ERROR == []),
identical to a genuine "no links" answer. Numeric UIDs also collide across databases
(PubMed 31322957 and nucleotide 31322957 are different records), so resolving is not
enough -- compare the printed label with the record you meant.

Usage: python check_source_ids.py <dbfrom> <uid> [<uid> ...]
       e.g. python check_source_ids.py gene 672 7157
Import: from check_source_ids import describe_ids, checked_elink
'''
import sys
import time
from Bio import Entrez

Entrez.email = 'your.email@example.com'
DELAY = 0.34


def describe_ids(dbfrom, ids):
    '''Return {uid: label}; raise ValueError naming every uid that does not resolve in dbfrom.'''
    ids = [str(i) for i in ids]
    labels, missing = {}, []
    for uid in ids:
        try:
            h = Entrez.esummary(db=dbfrom, id=uid)
            r = Entrez.read(h); h.close()
        except RuntimeError:  # "cannot get document summary"
            missing.append(uid)
            time.sleep(DELAY)
            continue
        # gene/pubmed/nuccore return a list of docsums; bioproject/sra/biosample return
        # {'DocumentSummarySet': {'DocumentSummary': [...]}}
        doc = r['DocumentSummarySet']['DocumentSummary'][0] if isinstance(r, dict) else r[0]
        parts = [str(doc.get(k)) for k in ('Caption', 'Name', 'Project_Acc', 'Title', 'Description') if doc.get(k)]
        labels[uid] = ' | '.join(parts)[:120] if parts else '<no label field>'
        time.sleep(DELAY)
    if missing:
        raise ValueError(f'UIDs not found in db={dbfrom!r}: {missing} -- wrong dbfrom, or an accession/symbol instead of a UID')
    return labels


def checked_elink(dbfrom, db, ids, **kwargs):
    '''ELink after confirming every source UID resolves in dbfrom; prints one label per UID.'''
    for uid, label in describe_ids(dbfrom, ids).items():
        print(f'  source {dbfrom}:{uid} = {label}')
    h = Entrez.elink(dbfrom=dbfrom, db=db, id=','.join(str(i) for i in ids), **kwargs)
    r = Entrez.read(h); h.close()
    return r


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    try:
        for uid, label in describe_ids(sys.argv[1], sys.argv[2:]).items():
            print(f'{sys.argv[1]}:{uid}\t{label}')
    except ValueError as e:
        sys.exit(f'ERROR: {e}')
