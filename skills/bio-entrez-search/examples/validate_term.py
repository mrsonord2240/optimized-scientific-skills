"""Local preflight for an Entrez query string; no network call is made."""


def validate_term(term, max_chars=2000, max_or_clauses=100):
    """Reject malformed or impractically large query strings before an API call."""
    if not isinstance(term, str) or not term.strip():
        raise ValueError('term must be a non-empty string')
    if len(term) > max_chars or term.upper().count(' OR ') > max_or_clauses:
        raise ValueError('term is too large; EPost IDs in batches instead')
    return term


if __name__ == '__main__':
    print(validate_term('BRCA1[Gene Name] AND Homo sapiens[Organism]'))
    for invalid in ('', 'BRCA1' * 1000, ' OR '.join(['BRCA1'] * 102)):
        try:
            validate_term(invalid)
        except ValueError as error:
            print(f'Rejected as expected: {error}')
        else:
            raise AssertionError(f'expected rejection: {invalid[:40]!r}')
