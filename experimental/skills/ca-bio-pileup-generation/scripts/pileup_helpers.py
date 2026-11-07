#!/usr/bin/env python3
"""Small, reusable pysam helpers for ``bio-pileup-generation``.

All coordinates accepted by these functions are 0-based.  The helpers cover
SNVs and ordinary aligner CIGARs; ``pileup_text`` deliberately refuses CIGAR
padding (``P``), whose mpileup rendering is not implemented here.
"""

from collections import Counter

import pysam


def allele_counts(bam_path, chrom, pos, **pileup_kw):
    """Return base counts (and ``DEL``) at 0-based *pos*.

    Reference skips and read-base ``N`` are excluded.  This is SNV-oriented:
    insertions or deletions following *pos* are not alleles at that position.
    Additional keyword arguments are passed to :meth:`AlignmentFile.pileup`.
    """
    counts = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for column in bam.pileup(chrom, pos, pos + 1, truncate=True, **pileup_kw):
            if column.pos != pos:
                continue
            for read in column.pileups:
                if read.is_refskip:
                    continue
                if read.is_del:
                    counts["DEL"] += 1
                    continue
                base = read.alignment.query_sequence[read.query_position].upper()
                if base != "N":
                    counts[base] += 1
    return dict(counts)


def allele_frequency(bam_path, chrom, pos, **pileup_kw):
    """Return non-deletion allele frequencies at 0-based *pos*."""
    counts = allele_counts(bam_path, chrom, pos, **pileup_kw)
    total = sum(count for base, count in counts.items() if base != "DEL")
    return {} if total == 0 else {
        base: count / total for base, count in counts.items() if base != "DEL"
    }


def find_variants(bam_path, ref_path, chrom, start, end, min_depth=10,
                  min_alt_freq=0.1, **pileup_kw):
    """Return SNVs against *ref_path* in 0-based, half-open ``[start, end)``.

    Reference-N sites, read-base N, and indels are excluded.  Base quality
    defaults to 20 unless the caller supplies ``min_base_quality``.
    """
    variants = []
    pileup_kw.setdefault("min_base_quality", 20)
    with (pysam.AlignmentFile(bam_path, "rb") as bam,
          pysam.FastaFile(ref_path) as reference):
        for column in bam.pileup(chrom, start, end, truncate=True, **pileup_kw):
            ref_base = reference.fetch(chrom, column.pos, column.pos + 1).upper()
            if ref_base == "N":
                continue
            alleles = Counter()
            for read in column.pileups:
                if read.is_refskip or read.is_del:
                    continue
                base = read.alignment.query_sequence[read.query_position].upper()
                if base != "N":
                    alleles[base] += 1
            depth = sum(alleles.values())
            if depth < min_depth:
                continue
            for base, count in alleles.items():
                if base != ref_base and count / depth >= min_alt_freq:
                    variants.append({"chrom": chrom, "pos": column.pos + 1,
                                     "ref": ref_base, "alt": base,
                                     "depth": depth, "alt_count": count,
                                     "freq": count / depth})
    return variants


def _indel_text(alignment, reference, chrom, pos):
    """Render the ``+N``/``-N`` marker following an aligned base."""
    ops = []
    for op, length in alignment.cigartuples:
        if ops and ops[-1][0] == op:
            ops[-1] = (op, ops[-1][1] + length)
        else:
            ops.append((op, length))
    ref_pos, query_pos = alignment.reference_start, 0
    for index, (op, length) in enumerate(ops):
        if op in (0, 2, 3, 7, 8):
            if ref_pos <= pos < ref_pos + length:
                break
            ref_pos += length
        if op in (0, 1, 4, 7, 8):
            query_pos += length
    if pos != ref_pos + length - 1:
        return ""
    if op in (0, 7, 8):
        query_pos += length
    remaining, text = ops[index + 1:], ""
    if remaining and remaining[0][0] == 1:
        length = remaining[0][1]
        text += f"+{length}" + alignment.query_sequence[query_pos:query_pos + length]
        remaining = remaining[1:]
    if remaining and remaining[0][0] == 2:
        length = remaining[0][1]
        text += f"-{length}" + reference.fetch(chrom, pos + 1, pos + 1 + length)
    return text.lower() if alignment.is_reverse else text.upper()


def pileup_text(bam_path, ref_path, chrom, start, end, **pileup_kw):
    """Yield six-column ``samtools mpileup -f``-style rows for ``[start, end)``.

    Pass ``compute_baq=False`` for ``-B``.  This helper supports ordinary
    aligner CIGARs but not padding (``P``); use ``samtools mpileup`` directly
    when padded alignments must be rendered.
    """
    options = {"stepper": "samtools", "truncate": True}
    options.update(pileup_kw)
    with (pysam.AlignmentFile(bam_path, "rb") as bam,
          pysam.FastaFile(ref_path) as reference):
        for column in bam.pileup(chrom, start, end, fastafile=reference, **options):
            ref_base = reference.fetch(chrom, column.pos, column.pos + 1)
            bases, qualities = [], []
            for read in column.pileups:
                alignment = read.alignment
                if any(op == 6 for op, _ in alignment.cigartuples):
                    raise ValueError("pileup_text does not support CIGAR padding (P); use samtools mpileup")
                text = ""
                if read.is_head:
                    text += "^" + chr(min(alignment.mapping_quality, 93) + 33)
                if read.is_refskip:
                    text += "<" if alignment.is_reverse else ">"
                elif read.is_del:
                    text += "*"
                else:
                    base = alignment.query_sequence[read.query_position]
                    if base.upper() == ref_base.upper():
                        text += "," if alignment.is_reverse else "."
                    else:
                        text += base.lower() if alignment.is_reverse else base.upper()
                text += _indel_text(alignment, reference, chrom, column.pos)
                if read.is_tail:
                    text += "$"
                bases.append(text)
                quality_pos = read.query_position_or_next
                if quality_pos is None or quality_pos >= len(alignment.query_qualities):
                    # samtools writes Q0 for a deletion after a read's last base.
                    qualities.append("!")
                else:
                    qualities.append(chr(min(alignment.query_qualities[quality_pos], 93) + 33))
            yield (f"{chrom}\t{column.pos + 1}\t{ref_base}\t{len(column.pileups)}\t"
                   f"{''.join(bases) or '*'}\t{''.join(qualities) or '*'}")
