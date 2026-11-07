'''Production-grade batch download: history server + disk checkpoint + WebEnv-expiry detection + jittered backoff.'''
# Verified: biopython 1.88, entrez direct 26.0 (2026-09-22) | Verify API if version differs
from Bio import Entrez
from urllib.error import HTTPError
import json
import random
import time
from pathlib import Path

Entrez.email = 'your.email@example.com'
# Entrez.api_key = 'your_api_key'


def truncate_to_offset(path, offset):
    '''Roll the output back to the last committed chunk boundary.

    The checkpoint records the output file's byte length at the moment the chunk
    was committed, so resuming is a byte-exact truncation -- no guessing where the
    previous process stopped.
    '''
    p = Path(path)
    if not p.exists():
        return
    with open(p, 'rb+') as f:
        f.truncate(min(offset, p.stat().st_size))


def truncate_to_last_newline(path):
    '''Fallback for checkpoints written before byte offsets were recorded.

    Drops the trailing partial line only. This is strictly weaker than
    truncate_to_offset: a crash after a chunk's records were written but before
    the checkpoint update leaves complete records past the committed boundary,
    and this cannot detect them.
    '''
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return
    with open(p, 'rb+') as f:
        f.seek(-min(8192, p.stat().st_size), 2)
        tail = f.read()
        last_nl = tail.rfind(b'\n')
        if last_nl >= 0:
            f.seek(-len(tail) + last_nl + 1, 2)
            f.truncate()


def checkpointed_download(db, term, out_path, ckpt_path, rettype='fasta',
                           batch_size=500, max_retries=5):
    delay = 0.1 if Entrez.api_key else 0.34
    ckpt = Path(ckpt_path)
    state = json.loads(ckpt.read_text()) if ckpt.exists() else {}
    start = state.get('start', 0)
    offset = state.get('offset')

    def refresh_session():
        h = Entrez.esearch(db=db, term=term, usehistory='y', retmax=0)
        s = Entrez.read(h); h.close()
        return s['WebEnv'], s['QueryKey'], int(s['Count'])

    webenv, query_key, total = refresh_session()
    print(f'{total:,} records matched; resuming at {start:,}')
    if start >= total:
        print('Already complete.')
        return

    if start == 0:
        Path(out_path).write_bytes(b'')
        offset = 0
    elif offset is None:
        print('  Checkpoint has no byte offset; falling back to last-newline truncation')
        truncate_to_last_newline(out_path)
        offset = Path(out_path).stat().st_size if Path(out_path).exists() else 0
    else:
        truncate_to_offset(out_path, offset)

    with open(out_path, 'ab') as out:
        while start < total:
            success = False
            for attempt in range(max_retries):
                try:
                    h = Entrez.efetch(db=db, rettype=rettype, retmode='text',
                                      retstart=start, retmax=batch_size,
                                      webenv=webenv, query_key=query_key)
                    body = h.read(); h.close()
                    raw = body if isinstance(body, bytes) else body.encode('utf-8')
                    text = raw.decode('utf-8', errors='replace')
                    if not text.strip():
                        raise RuntimeError('Empty body')
                    if '<ERROR>' in text[:500]:
                        raise RuntimeError(f'WebEnv expired or server error: {text[:200]}')
                    out.write(raw); out.flush()
                    success = True
                    break
                except HTTPError as e:
                    backoff = min(120, (2 ** attempt) + random.uniform(0, 1))
                    if e.code == 429:
                        print(f'  HTTP 429; backing off {backoff:.1f}s')
                        time.sleep(backoff)
                    elif e.code in (500, 502, 503, 504):
                        print(f'  HTTP {e.code}; backing off {backoff:.1f}s')
                        time.sleep(backoff)
                    else:
                        raise
                except RuntimeError as e:
                    print(f'  {e}; refreshing WebEnv')
                    webenv, query_key, total = refresh_session()
                    time.sleep(1)

            if not success:
                raise RuntimeError(f'Failed after {max_retries} retries at start={start}')

            start += batch_size
            # Byte offset of the committed output, written together with the cursor.
            # Both are flushed to disk before the checkpoint lands, so a crash can
            # only ever leave *extra* bytes past `offset` -- never a checkpoint
            # pointing past what was actually written.
            offset = out.tell()
            ckpt.write_text(json.dumps({'start': start, 'total': total, 'offset': offset}))
            time.sleep(delay)
            if (start // batch_size) % 10 == 0 or start >= total:
                print(f'  {min(start, total):,}/{total:,}')

    ckpt.unlink(missing_ok=True)
    print(f'Done: {total:,} records -> {out_path}')


if __name__ == '__main__':
    checkpointed_download(
        db='nucleotide',
        # Gene-symbol field tags need the official gene symbol, not a descriptive word
        # (checked live 2026-09-19: 'hemoglobin[Gene Name]' returns Count=0).
        term='BRCA1[GENE] AND Homo sapiens[ORGN] AND biomol_mrna[PROP] AND srcdb_refseq[PROP]',
        out_path='brca1_mrna.fasta',
        ckpt_path='brca1_mrna.ckpt.json',
    )