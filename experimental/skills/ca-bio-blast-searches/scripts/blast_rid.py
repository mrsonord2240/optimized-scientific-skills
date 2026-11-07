'''Submit a remote BLAST search, capture its RID, poll independently, fetch the XML (NCBI BLAST URL API, stdlib only).

Inputs : FASTA file with a defline (>id), program, database; optional hitlist size / expect / megablast flags.
Usage  : python blast_rid.py run    --query q.fa --program blastn --db refseq_select_rna --out hits.xml
         python blast_rid.py submit --query q.fa --program blastn --db refseq_select_rna      # prints RID and RTOE
         python blast_rid.py status RID
         python blast_rid.py fetch  RID --out hits.xml
Checked : NCBI BLAST URL API 2026-09-21; Biopython 1.88 NCBIXML.read() parses the fetched XML.
Polite  : waits RTOE before the first poll, then polls no more than once per 60 s (NCBI's stated limit).
'''
import argparse
import http.client
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'
TOOL = {'TOOL': 'bioskills-blast-searches'}  # add 'EMAIL': you@example.org so NCBI can contact you


def _call(params, post=False):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=data) if post else urllib.request.Request(f'{URL}?{data.decode()}')
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.read().decode('utf-8', 'replace')
    except (urllib.error.URLError, http.client.HTTPException, OSError) as e:  # network error or NCBI maintenance: retry later, the RID stays valid
        sys.exit(f'NCBI request failed: {e}')


def submit(query_fasta, program, db, hitlist_size=500, expect=None, megablast=False, extra=None):
    if not query_fasta.lstrip().startswith('>'):
        sys.exit('Query needs a FASTA defline (>id): anonymous queries run ~12x slower (see SKILL.md Failure Modes)')
    params = {'CMD': 'Put', 'PROGRAM': program, 'DATABASE': db, 'QUERY': query_fasta,
              'HITLIST_SIZE': hitlist_size, **TOOL}
    if expect is not None:
        params['EXPECT'] = expect
    if megablast:
        params['MEGABLAST'] = 'on'
    params.update(extra or {})
    html = _call(params, post=True)
    rid = re.search(r'RID = (\S+)', html)
    rtoe = re.search(r'RTOE = (\d+)', html)
    if not rid:
        msg = re.search(r'Message ID#(\d+) Error: ([^<\n]+)', html)
        sys.exit(f'Submit rejected: {msg.group(0) if msg else html[:300]}')
    return rid.group(1), int(rtoe.group(1)) if rtoe else 15


def status(rid):
    '''Return (Status, ThereAreHits): Status is WAITING, READY, FAILED or UNKNOWN (UNKNOWN = expired/invalid RID).'''
    html = _call({'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid})
    st = re.search(r'Status=(\w+)', html)
    return (st.group(1) if st else 'UNKNOWN'), 'ThereAreHits=yes' in html


def fetch(rid, out_path):
    xml = _call({'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'XML'})
    if '<BlastOutput>' not in xml[:2000] and '<BlastOutput ' not in xml[:2000]:
        sys.exit(f'No BLAST XML returned for {rid} (not READY yet, RID expired after 24-36 h, or search failed; check status first): {xml[:200]!r}')
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(xml)
    return out_path


def wait(rid, rtoe, interval=60, timeout=3600):
    '''Poll until READY. Permissive searches (PAM30, word_size=2, expect=1000, huge hitlist) and busy queues can take
    many minutes (one run took 19 min): a long WAITING is not a hang. Raises on FAILED/UNKNOWN or timeout.'''
    t0 = time.time()
    time.sleep(rtoe)
    while True:
        st, hits = status(rid)
        print(f'[{time.time() - t0:6.0f}s] {rid} {st}', flush=True)
        if st == 'READY':
            return hits
        if st in ('FAILED', 'UNKNOWN'):
            sys.exit(f'{rid} {st}: FAILED = NCBI-side error (resubmit, smaller query); UNKNOWN = expired/invalid RID')
        if time.time() - t0 > timeout:
            sys.exit(f'Timed out after {timeout}s; {rid} is still valid, resume with: status {rid} / fetch {rid}')
        time.sleep(interval)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    for name in ('run', 'submit'):
        p = sub.add_parser(name)
        p.add_argument('--query', required=True, help='FASTA file (with defline)')
        p.add_argument('--program', default='blastn')
        p.add_argument('--db', default='refseq_select_rna')
        p.add_argument('--hitlist', type=int, default=500)
        p.add_argument('--expect', type=float)
        p.add_argument('--megablast', action='store_true', help='blastn + MEGABLAST=on')
        if name == 'run':
            p.add_argument('--out', required=True)
            p.add_argument('--interval', type=int, default=60)
            p.add_argument('--timeout', type=int, default=3600)
    p = sub.add_parser('status'); p.add_argument('rid')
    p = sub.add_parser('fetch'); p.add_argument('rid'); p.add_argument('--out', required=True)
    a = ap.parse_args()

    if a.cmd in ('run', 'submit'):
        with open(a.query, encoding='utf-8') as f:
            rid, rtoe = submit(f.read(), a.program, a.db, a.hitlist, a.expect, a.megablast)
        print(f'RID={rid} RTOE={rtoe}s  https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}', flush=True)
        if a.cmd == 'run':
            hits = wait(rid, rtoe, a.interval, a.timeout)
            print('ThereAreHits=' + ('yes' if hits else 'no'))
            print('saved', fetch(rid, a.out))
    elif a.cmd == 'status':
        print(*status(a.rid))
    else:
        print('saved', fetch(a.rid, a.out))


if __name__ == '__main__':
    main()
