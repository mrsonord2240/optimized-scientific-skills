"""Safe pybiomart 0.2.0 XML querying for BioMart ID-list filters.

Inputs: a pybiomart Dataset, attribute names, and filters. Usage:
    from biomart_query import query_raw
    frame = query_raw(dataset, attributes=[...], filters={"ensembl_gene_id": ids})

Run ``python scripts/biomart_query.py --self-test`` for an offline validation
of batching and input-validation rules. The helper supports one large
list-valued filter at a time; pre-filter independent large lists upstream.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from io import StringIO
from xml.etree import ElementTree

import pandas as pd
from requests.exceptions import HTTPError


DEFAULT_ID_LIST_CHUNK_SIZE = 500


class BioMartInputError(ValueError):
    """The request is invalid locally and was not sent to BioMart."""


class BioMartResponseError(RuntimeError):
    """BioMart responded, but not with a parseable TSV result."""


class BioMartOutageError(BioMartResponseError):
    """BioMart returned an HTML/status response that can be retried later."""


def _as_values(value: object) -> list[object] | None:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        return None
    return list(value)


def _query_once(ds, attributes: Sequence[str], filters: Mapping[str, object]):
    root = ElementTree.Element("Query")
    root.set("virtualSchemaName", "default")
    root.set("formatter", "TSV")
    root.set("header", "1")
    root.set("uniqueRows", "1")
    root.set("datasetConfigVersion", "0.6")
    dataset_el = ElementTree.SubElement(root, "Dataset")
    dataset_el.set("name", ds.name)
    dataset_el.set("interface", "default")
    for name, value in filters.items():
        values = _as_values(value)
        filter_el = ElementTree.SubElement(dataset_el, "Filter")
        filter_el.set("name", name)
        filter_el.set("value", ",".join(map(str, values)) if values is not None else str(value))
    for name in attributes:
        attribute_el = ElementTree.SubElement(dataset_el, "Attribute")
        attribute_el.set("name", name)

    body = ds.get(query=ElementTree.tostring(root)).text.strip()
    body_lower = body.lower()
    if "query error" in body_lower:
        raise BioMartResponseError(f"BioMart rejected the query: {body}")
    if not body:
        raise BioMartResponseError(
            "BioMart returned an empty body; retry once, then check filters and service status."
        )
    if (body_lower.startswith("<html") or body_lower.startswith("<!doctype html")
            or "status.ensembl.org" in body_lower or "service unavailable" in body_lower):
        raise BioMartOutageError(
            "BioMart returned an HTML/status outage page instead of TSV; retry with backoff."
        )
    try:
        return pd.read_csv(StringIO(body), sep="\t")
    except Exception as exc:
        raise BioMartResponseError(
            "BioMart returned non-TSV text; inspect the response before retrying."
        ) from exc


def query_raw(ds, attributes: Sequence[str], filters: Mapping[str, object],
              *, id_list_chunk_size: int = DEFAULT_ID_LIST_CHUNK_SIZE):
    """Query BioMart despite pybiomart's ID-list validation gap.

    A non-empty ID list longer than ``id_list_chunk_size`` is split into
    independent requests and concatenated. If an instance has a shorter GET
    URL limit, a 414 response halves only the failed batch and retries it.
    """
    if not attributes:
        raise BioMartInputError("attributes must name at least one BioMart field")
    if not filters:
        raise BioMartInputError("filters must contain at least one BioMart constraint")
    if id_list_chunk_size < 1:
        raise BioMartInputError("id_list_chunk_size must be a positive integer")

    list_filters = []
    for name, value in filters.items():
        values = _as_values(value)
        if values is not None:
            if not values:
                raise BioMartInputError(f"filter '{name}' is empty; no request was sent")
            if len(values) > id_list_chunk_size:
                list_filters.append((name, values))
    if len(list_filters) > 1:
        names = ", ".join(name for name, _ in list_filters)
        raise BioMartInputError(
            f"only one list-valued filter may exceed {id_list_chunk_size} items ({names}); "
            "pre-filter one list before querying BioMart"
        )
    if not list_filters:
        return _query_once(ds, attributes, filters)

    filter_name, values = list_filters[0]
    pending_batches = [values[start:start + id_list_chunk_size]
                       for start in range(0, len(values), id_list_chunk_size)]
    frames = []
    while pending_batches:
        batch_values = pending_batches.pop(0)
        batch_filters = dict(filters)
        batch_filters[filter_name] = batch_values
        try:
            frames.append(_query_once(ds, attributes, batch_filters))
        except HTTPError as exc:
            status = getattr(exc.response, "status_code", None)
            if status != 414 or len(batch_values) == 1:
                raise
            midpoint = len(batch_values) // 2
            pending_batches[:0] = [batch_values[:midpoint], batch_values[midpoint:]]
    return pd.concat(frames, ignore_index=True)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Validate BioMart batching rules without a network call.")
    parser.add_argument("--self-test", action="store_true", help="run offline helper checks")
    args = parser.parse_args()
    if not args.self_test:
        parser.error("use --self-test, or import query_raw from this module")
    assert _as_values("ENSG000001") is None
    assert _as_values(["A", "B"]) == ["A", "B"]
    try:
        query_raw(None, ["ensembl_gene_id"], {"ensembl_gene_id": []})
    except BioMartInputError:
        pass
    else:
        raise AssertionError("empty ID lists must fail before any network call")
    print("offline BioMart helper checks passed")
