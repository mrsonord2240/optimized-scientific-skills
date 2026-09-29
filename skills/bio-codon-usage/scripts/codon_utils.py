"""Shared, table-aware validation and CAI helpers for codon-usage analyses."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Literal

from Bio.Data import CodonTable
from Bio.Seq import Seq


ValidationPolicy = Literal["strict", "permissive"]


class CDSValidationError(ValueError):
    """Raised when a sequence cannot satisfy the requested CDS contract."""


class CAIScoringError(ValueError):
    """Raised when a validated query has no codons included in a CAI score."""


@dataclass(frozen=True)
class DiscardedSegment:
    offset: int
    text: str
    reason: str


@dataclass(frozen=True)
class CDSValidationResult:
    sequence: str
    codons: tuple[str, ...]
    table_id: int
    policy: ValidationPolicy
    starts_with_allowed_codon: bool
    has_terminal_stop: bool
    internal_stop_indices: tuple[int, ...]
    discarded: tuple[DiscardedSegment, ...]

    def discard_summary(self) -> str:
        if not self.discarded:
            return "none"
        return "; ".join(
            f"offset {item.offset}: {item.text!r} ({item.reason})"
            for item in self.discarded
        )


@dataclass(frozen=True)
class CAIScoreResult:
    score: float
    included_codons: int
    excluded_codons: tuple[tuple[int, str, str], ...]
    validation: CDSValidationResult


def _table(table_id: int):
    try:
        return CodonTable.unambiguous_dna_by_id[table_id]
    except KeyError as exc:
        raise CDSValidationError(f"unknown NCBI genetic-code table: {table_id}") from exc


def validate_cds(
    sequence,
    *,
    table_id: int = 1,
    policy: ValidationPolicy = "strict",
    require_start: bool = True,
    require_terminal_stop: bool = True,
    allow_internal_stops: bool = False,
) -> CDSValidationResult:
    """Validate a DNA CDS and return the exact sequence admitted to analysis.

    Strict mode rejects non-ACGT triplets and incomplete trailing bases.
    Permissive mode discards whole invalid triplets and a trailing 1-2 base
    remainder, recording every discard. Both modes apply the requested start,
    stop, and internal-stop checks using the selected genetic-code table.
    """

    if policy not in {"strict", "permissive"}:
        raise ValueError("policy must be 'strict' or 'permissive'")

    raw = str(sequence).upper()
    if not raw:
        raise CDSValidationError("CDS is empty; provide at least one complete codon")
    if any(character.isspace() for character in raw):
        raise CDSValidationError("CDS contains whitespace; normalize FASTA input before validation")

    remainder = len(raw) % 3
    discarded: list[DiscardedSegment] = []
    complete = raw
    if remainder:
        offset = len(raw) - remainder
        trailing = raw[offset:]
        if policy == "strict":
            raise CDSValidationError(
                f"CDS length {len(raw)} is not divisible by three; trailing bases {trailing!r} at offset {offset}"
            )
        discarded.append(DiscardedSegment(offset, trailing, "incomplete trailing codon"))
        complete = raw[:offset]

    retained: list[str] = []
    for offset in range(0, len(complete), 3):
        codon = complete[offset : offset + 3]
        illegal = sorted(set(codon) - set("ACGT"))
        if illegal:
            if policy == "strict":
                raise CDSValidationError(
                    f"invalid DNA codon {codon!r} at offset {offset}; illegal symbols: {''.join(illegal)}"
                )
            discarded.append(DiscardedSegment(offset, codon, "codon contains non-ACGT symbols"))
            continue
        retained.append(codon)

    if not retained:
        raise CDSValidationError("no complete unambiguous codons remain after validation")

    discarded.sort(key=lambda item: item.offset)

    table = _table(table_id)
    stop_indices = tuple(index for index, codon in enumerate(retained) if codon in table.stop_codons)
    internal_stops = tuple(index for index in stop_indices if index != len(retained) - 1)
    starts_correctly = retained[0] in table.start_codons
    terminal_stop = retained[-1] in table.stop_codons

    if require_start and not starts_correctly:
        raise CDSValidationError(
            f"first codon {retained[0]} is not a start codon in genetic-code table {table_id}"
        )
    if require_terminal_stop and not terminal_stop:
        raise CDSValidationError(
            f"last codon {retained[-1]} is not a stop codon in genetic-code table {table_id}"
        )
    if internal_stops and not allow_internal_stops:
        raise CDSValidationError(
            f"internal stop codon(s) at codon indices {list(internal_stops)} under genetic-code table {table_id}"
        )

    return CDSValidationResult(
        sequence="".join(retained),
        codons=tuple(retained),
        table_id=table_id,
        policy=policy,
        starts_with_allowed_codon=starts_correctly,
        has_terminal_stop=terminal_stop,
        internal_stop_indices=internal_stops,
        discarded=tuple(discarded),
    )


def count_codons(sequence, **validation_options):
    """Return validated codon counts and the validation/discard report."""

    validation = validate_cds(sequence, **validation_options)
    return Counter(validation.codons), validation


def codon_frequencies(sequence, **validation_options):
    """Return validated codon frequencies and the validation/discard report."""

    counts, validation = count_codons(sequence, **validation_options)
    total = sum(counts.values())
    return {codon: count / total for codon, count in counts.items()}, validation


def _assert_cai_table(cai, table_id: int) -> None:
    selected = _table(table_id)
    actual = getattr(cai, "_table", None)
    if actual is None:
        raise ValueError("CAI object does not expose its genetic-code table; cannot verify table consistency")
    if actual.forward_table != selected.forward_table or set(actual.stop_codons) != set(selected.stop_codons):
        raise ValueError(f"CAI index was not constructed with genetic-code table {table_id}")


def calculate_cai_guarded(
    cai,
    sequence,
    *,
    table_id: int = 1,
    policy: ValidationPolicy = "strict",
    require_start: bool = True,
    require_terminal_stop: bool = True,
) -> CAIScoreResult:
    """Validate and score a CDS, guarding Biopython's zero denominator.

    Biopython 1.85 always excludes ATG and TGG. Standard stops absent from the
    index are also skipped; any codon present in the mapping, including a stop,
    contributes. This wrapper reports those exclusions before calling the live
    implementation.
    """

    _assert_cai_table(cai, table_id)
    validation = validate_cds(
        sequence,
        table_id=table_id,
        policy=policy,
        require_start=require_start,
        require_terminal_stop=require_terminal_stop,
    )
    excluded: list[tuple[int, str, str]] = []
    included = 0
    for index, codon in enumerate(validation.codons):
        if codon in {"ATG", "TGG"}:
            excluded.append((index, codon, "Biopython 1.85 hard-coded exclusion"))
        elif codon not in cai and codon in {"TAA", "TAG", "TGA"}:
            excluded.append((index, codon, "standard stop absent from CAI mapping"))
        else:
            included += 1
    if included == 0:
        excluded_text = ", ".join(f"{index}:{codon}" for index, codon, _ in excluded)
        raise CAIScoringError(
            "CAI is undefined because the validated query has no included codons"
            + (f" (excluded: {excluded_text})" if excluded_text else "")
        )
    score = float(cai.calculate(Seq(validation.sequence)))
    return CAIScoreResult(score, included, tuple(excluded), validation)


def optimize_dna_table_aware(
    cai,
    sequence,
    *,
    table_id: int = 1,
    strict: bool = True,
    policy: ValidationPolicy = "strict",
):
    """Optimize DNA through an explicitly table-aware protein intermediate."""

    _assert_cai_table(cai, table_id)
    validation = validate_cds(sequence, table_id=table_id, policy=policy)
    protein = Seq(validation.sequence).translate(table=table_id)
    optimized = cai.optimize(protein, seq_type="protein", strict=strict)
    optimized_validation = validate_cds(optimized, table_id=table_id, policy="strict")
    optimized_protein = Seq(optimized_validation.sequence).translate(table=table_id)
    if optimized_protein != protein:
        raise RuntimeError(
            f"table-aware optimization changed the protein under genetic-code table {table_id}"
        )
    return optimized, validation, optimized_validation
