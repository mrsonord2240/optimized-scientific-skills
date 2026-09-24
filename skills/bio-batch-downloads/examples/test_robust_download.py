'''Regression tests for checkpointed_download: resume correctness and the error guards.

No network. The fake Entrez stands in for the history server and serves records
from a locally synthesized set, sliced at record boundaries exactly as NCBI does.

Run: python -m pytest examples/test_robust_download.py -q
'''
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import robust_download as rd


# --- fake NCBI history server -------------------------------------------------

class _Handle:
    def __init__(self, payload):
        self._payload = payload

    def read(self):
        return self._payload

    def close(self):
        pass


class FakeEntrez:
    '''Serves a fixed record set; efetch slices it by retstart/retmax.'''
    api_key = 'test-key'  # shortens the inter-chunk sleep

    def __init__(self, records):
        self.records = records
        self.efetch_calls = []

    def esearch(self, **kwargs):
        return _Handle({'WebEnv': 'WEBENV', 'QueryKey': '1', 'Count': str(len(self.records))})

    def efetch(self, **kwargs):
        start, n = kwargs['retstart'], kwargs['retmax']
        self.efetch_calls.append(start)
        return _Handle(b''.join(self.records[start:start + n]))

    def read(self, handle):
        return handle.read()


def make_records(n, seq_len=120):
    '''n synthetic FASTA records, each terminated by a newline.

    Records are separated by a blank line, matching NCBI's actual EFetch FASTA
    layout (`...TCCA\\n\\n>NM_...`). The byte offset arithmetic only holds if the
    separator is reproduced, so the fixture must not normalize it away.
    '''
    recs = [
        f'>NM_{i:06d}.1 synthetic record {i}\n'.encode() + (b'ACGT' * (seq_len // 4)) + b'\n'
        for i in range(n)
    ]
    return [r + b'\n' for r in recs[:-1]] + recs[-1:]


@pytest.fixture
def fake(monkeypatch):
    records = make_records(368)
    f = FakeEntrez(records)
    monkeypatch.setattr(rd, 'Entrez', f)
    return f


# --- the resume defect --------------------------------------------------------

def test_resume_after_midchunk_crash_is_byte_exact(tmp_path, fake):
    '''A crash that wrote past the committed boundary must not duplicate records.

    Reproduces the real failure: 150 complete records on disk plus the torn start
    of #151, while the checkpoint only vouches for 100. Truncating to the last
    newline leaves the dangling header and the extra 50 records, so the resumed
    file ends up longer than the reference.
    '''
    records = fake.records
    out = tmp_path / 'mid.fasta'
    ckpt = tmp_path / 'mid.ckpt.json'

    committed = 100
    offset = len(b''.join(records[:committed]))
    torn = b''.join(records[:150]) + records[150].split(b'\n')[0] + b'\n' + b'ACGTA'
    out.write_bytes(torn)
    ckpt.write_text(json.dumps({'start': committed, 'total': len(records), 'offset': offset}))

    rd.checkpointed_download('nucleotide', 'TERM', str(out), str(ckpt), batch_size=100)

    assert out.read_bytes() == b''.join(records), 'resumed file is not byte-identical to the reference'
    assert not ckpt.exists(), 'checkpoint should be removed on success'


def test_resume_rolls_back_to_committed_offset(tmp_path, fake):
    '''Bytes past the committed offset are discarded even when they look complete.'''
    records = fake.records
    out = tmp_path / 'rollback.fasta'
    ckpt = tmp_path / 'rollback.ckpt.json'

    committed = 200
    offset = len(b''.join(records[:committed]))
    # 40 *complete* records past the boundary -- invisible to a newline heuristic
    out.write_bytes(b''.join(records[:240]))
    ckpt.write_text(json.dumps({'start': committed, 'total': len(records), 'offset': offset}))

    rd.checkpointed_download('nucleotide', 'TERM', str(out), str(ckpt), batch_size=100)

    assert out.read_bytes() == b''.join(records)


def test_legacy_checkpoint_without_offset_still_completes(tmp_path, fake):
    '''Old-format checkpoints fall back to last-newline truncation and finish.'''
    records = fake.records
    out = tmp_path / 'legacy.fasta'
    ckpt = tmp_path / 'legacy.ckpt.json'

    out.write_bytes(b''.join(records[:100]))  # clean boundary, no torn tail
    ckpt.write_text(json.dumps({'start': 100, 'total': len(records)}))

    rd.checkpointed_download('nucleotide', 'TERM', str(out), str(ckpt), batch_size=100)

    assert out.read_bytes() == b''.join(records)


# --- guards -------------------------------------------------------------------

def test_zero_count_writes_no_file(tmp_path, monkeypatch):
    f = FakeEntrez([])
    monkeypatch.setattr(rd, 'Entrez', f)
    out = tmp_path / 'none.fasta'

    rd.checkpointed_download('nucleotide', 'TERM', str(out), str(tmp_path / 'none.ckpt.json'))

    assert not out.exists(), 'a zero-result query must not leave an empty file behind'
    assert f.efetch_calls == [], 'a zero-result query must not call efetch'


def test_str_body_is_decoded_not_crashed(tmp_path, monkeypatch):
    '''A str body (some Entrez paths return text, not bytes) is handled.'''
    records = make_records(20)

    class StrEntrez(FakeEntrez):
        def efetch(self, **kwargs):
            start, n = kwargs['retstart'], kwargs['retmax']
            self.efetch_calls.append(start)
            return _Handle(b''.join(self.records[start:start + n]).decode())

    f = StrEntrez(records)
    monkeypatch.setattr(rd, 'Entrez', f)
    out = tmp_path / 'str.fasta'

    rd.checkpointed_download('nucleotide', 'TERM', str(out), str(tmp_path / 'str.ckpt.json'))

    assert out.read_bytes() == b''.join(records)


def test_error_body_raises_after_retries(tmp_path, monkeypatch):
    '''HTTP 200 with an <ERROR> body (expired WebEnv) must not be written to disk.'''
    class ErrorEntrez(FakeEntrez):
        def efetch(self, **kwargs):
            self.efetch_calls.append(kwargs['retstart'])
            return _Handle(b'<ERROR>WebEnv not found</ERROR>')

    f = ErrorEntrez(make_records(50))
    monkeypatch.setattr(rd, 'Entrez', f)
    out = tmp_path / 'err.fasta'

    with pytest.raises(RuntimeError, match='retries'):
        rd.checkpointed_download('nucleotide', 'TERM', str(out),
                                 str(tmp_path / 'err.ckpt.json'), batch_size=10, max_retries=2)

    assert out.read_bytes() == b'', 'no error payload may reach the output file'
    assert not (tmp_path / 'err.ckpt.json').exists()


def test_completed_checkpoint_is_a_noop(tmp_path, fake):
    '''Resuming a finished job must not re-download or rewrite the output.'''
    out = tmp_path / 'done.fasta'
    out.write_bytes(b''.join(fake.records))
    before = out.read_bytes()
    ckpt = tmp_path / 'done.ckpt.json'
    ckpt.write_text(json.dumps({'start': len(fake.records), 'total': len(fake.records)}))

    rd.checkpointed_download('nucleotide', 'TERM', str(out), str(ckpt))

    assert out.read_bytes() == before
    assert fake.efetch_calls == []