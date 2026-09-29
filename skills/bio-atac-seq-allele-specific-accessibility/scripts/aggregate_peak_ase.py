#!/usr/bin/env python3
"""Aggregate phased GATK ASE counts within ATAC peaks.

REF is a locus-specific label and must never be summed across loci. This
helper therefore requires an explicit map from each variant to haplotype 1 and
keeps phase blocks separate when a peak crosses a block boundary.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pybedtools
from scipy import stats


ASE_COLUMNS = {
    "contig", "position", "variantID", "refAllele", "altAllele",
    "refCount", "altCount", "totalCount",
}
MAP_COLUMNS = {
    "contig", "position", "refAllele", "altAllele", "phaseSet",
    "haplotype1Allele",
}
OUTPUT_COLUMNS = [
    "peak", "phase_set", "haplotype1_count", "haplotype2_count",
    "total_count", "haplotype1_frac", "snp_count", "p_value", "adj_p",
    "status",
]
PEAK_ID_COLUMNS = ["peak_chrom", "peak_start", "peak_end"]


class InputContractError(ValueError):
    """An input is readable but violates the documented analysis contract."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Pool phased heterozygous SNP counts within ATAC peaks, oriented "
            "to haplotype 1 and separated by phase set."
        )
    )
    parser.add_argument("--ase-counts", required=True, help="GATK ASEReadCounter TABLE")
    parser.add_argument("--peaks", required=True, help="Consensus peaks BED3 or wider")
    parser.add_argument(
        "--haplotype-map", required=True,
        help=(
            "TSV with contig, position, refAllele, altAllele, phaseSet, and "
            "haplotype1Allele (REF or ALT)"
        ),
    )
    parser.add_argument("--output", required=True, help="All evaluated peak/phase-set TSV")
    return parser.parse_args()


def _require_columns(frame: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise InputContractError(f"{label} is missing required columns: {', '.join(missing)}")


def _integer_column(frame: pd.DataFrame, column: str, label: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.isna().any() or not np.equal(values, np.floor(values)).all():
        raise InputContractError(f"{label} column {column!r} must contain integers")
    return values.astype("int64")


def read_ase(path: str) -> pd.DataFrame:
    try:
        ase = pd.read_csv(path, sep="\t", dtype={"contig": str})
    except Exception as exc:
        raise InputContractError(f"could not read ASE table {path!r}: {exc}") from exc
    _require_columns(ase, ASE_COLUMNS, "ASE table")
    if ase.empty:
        return ase
    for column in ("position", "refCount", "altCount", "totalCount"):
        ase[column] = _integer_column(ase, column, "ASE table")
    if (ase[["position", "refCount", "altCount", "totalCount"]] < 0).any().any():
        raise InputContractError("ASE positions and counts must be non-negative")
    if (ase["refCount"] + ase["altCount"] > ase["totalCount"]).any():
        raise InputContractError("ASE refCount + altCount cannot exceed totalCount")
    keys = ["contig", "position", "refAllele", "altAllele"]
    if ase.duplicated(keys).any():
        raise InputContractError("ASE table contains duplicate variant keys")
    return ase


def read_haplotype_map(path: str) -> pd.DataFrame:
    try:
        phase = pd.read_csv(path, sep="\t", dtype={"contig": str, "phaseSet": str})
    except Exception as exc:
        raise InputContractError(f"could not read haplotype map {path!r}: {exc}") from exc
    _require_columns(phase, MAP_COLUMNS, "haplotype map")
    if phase.empty:
        raise InputContractError("haplotype map has no variants")
    phase["position"] = _integer_column(phase, "position", "haplotype map")
    if (phase["position"] < 1).any():
        raise InputContractError("haplotype-map positions must be one-based positive integers")
    if phase["phaseSet"].isna().any() or (phase["phaseSet"].str.strip() == "").any():
        raise InputContractError("every haplotype-map row needs a non-empty phaseSet")
    allowed = {"REF", "ALT"}
    observed = set(phase["haplotype1Allele"].dropna())
    if not observed.issubset(allowed) or phase["haplotype1Allele"].isna().any():
        raise InputContractError("haplotype1Allele must be REF or ALT on every row")
    keys = ["contig", "position", "refAllele", "altAllele"]
    if phase.duplicated(keys).any():
        raise InputContractError("haplotype map contains duplicate variant keys")
    return phase


def read_peaks(path: str) -> pd.DataFrame:
    rows: list[list[str]] = []
    with open(path, encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if not rows and len(fields) >= 3 and fields[0].lower() in {"chrom", "chr", "contig"}:
                continue
            if len(fields) < 3:
                raise InputContractError(f"peaks line {line_number} has fewer than three columns")
            rows.append(fields)
    peaks = pd.DataFrame(rows)
    if peaks.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "peak"])
    starts = pd.to_numeric(peaks[1], errors="coerce")
    ends = pd.to_numeric(peaks[2], errors="coerce")
    if starts.isna().any() or ends.isna().any():
        raise InputContractError("peak start and end columns must contain integers")
    if not np.equal(starts, np.floor(starts)).all() or not np.equal(ends, np.floor(ends)).all():
        raise InputContractError("peak start and end columns must contain integers")
    starts = starts.astype("int64")
    ends = ends.astype("int64")
    if (starts < 0).any() or (ends <= starts).any():
        raise InputContractError("peaks require 0 <= start < end")
    names = []
    for index, row in peaks.iterrows():
        supplied = str(row[3]).strip() if len(row) > 3 and pd.notna(row[3]) else ""
        names.append(supplied or f"{row[0]}_{starts[index]}_{ends[index]}")
    result = pd.DataFrame(
        {"chrom": peaks[0].astype(str), "start": starts, "end": ends, "peak": names}
    )
    if result.duplicated(["chrom", "start", "end"]).any():
        raise InputContractError("peaks contain duplicate genomic intervals")
    return result


def empty_output(path: str) -> None:
    pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(path, sep="\t", index=False)


def map_variants_to_peaks(ase: pd.DataFrame, peaks: pd.DataFrame) -> pd.DataFrame:
    if ase.empty or peaks.empty:
        return pd.DataFrame(columns=list(ase.columns) + PEAK_ID_COLUMNS + ["peak"])
    snps = pd.DataFrame({
        "chrom": ase["contig"], "start": ase["position"] - 1,
        "end": ase["position"], "variant_row": np.arange(len(ase)),
    })
    intersected = pybedtools.BedTool.from_dataframe(snps).intersect(
        pybedtools.BedTool.from_dataframe(peaks), wa=True, wb=True
    )
    if intersected.count() == 0:
        return pd.DataFrame(columns=list(ase.columns) + PEAK_ID_COLUMNS + ["peak"])
    mapped = intersected.to_dataframe(
        disable_auto_names=True, header=None,
        names=[
            "snp_chrom", "snp_start", "snp_end", "variant_row",
            "peak_chrom", "peak_start", "peak_end", "peak",
        ],
    )
    mapped["variant_row"] = pd.to_numeric(mapped["variant_row"], errors="raise").astype(int)
    return ase.iloc[mapped["variant_row"].to_numpy()].reset_index(drop=True).assign(
        peak_chrom=mapped["peak_chrom"].astype(str).to_numpy(),
        peak_start=pd.to_numeric(mapped["peak_start"], errors="raise").astype("int64").to_numpy(),
        peak_end=pd.to_numeric(mapped["peak_end"], errors="raise").astype("int64").to_numpy(),
        peak=mapped["peak"].astype(str).to_numpy(),
    )


def evaluate_groups(mapped: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    group_columns = PEAK_ID_COLUMNS + ["phaseSet"]
    for (_, _, _, phase_set), group in mapped.groupby(group_columns, sort=True):
        # BED4 names are optional display metadata. Genomic coordinates define
        # statistical peak identity, so repeated labels cannot pool loci.
        peak = str(group["peak"].iloc[0])
        hap1 = int(group["haplotype1Count"].sum())
        hap2 = int(group["haplotype2Count"].sum())
        total = hap1 + hap2
        snp_count = len(group)
        p_value = stats.binomtest(hap1, total, p=0.5).pvalue if snp_count >= 2 and total > 0 else np.nan
        rows.append({
            "peak": peak, "phase_set": phase_set,
            "haplotype1_count": hap1, "haplotype2_count": hap2,
            "total_count": total,
            "haplotype1_frac": hap1 / total if total else np.nan,
            "snp_count": snp_count, "p_value": p_value,
        })
    result = pd.DataFrame(rows)
    result["adj_p"] = np.nan
    tested = result["p_value"].notna()
    if tested.any():
        result.loc[tested, "adj_p"] = stats.false_discovery_control(
            result.loc[tested, "p_value"].to_numpy()
        )
    result["status"] = "non_significant"
    result.loc[result["snp_count"] < 2, "status"] = "underpowered_snp_count"
    result.loc[
        tested & (result["adj_p"] < 0.05) & (abs(result["haplotype1_frac"] - 0.5) < 0.2), "status"
    ] = "below_effect_threshold"
    result.loc[
        tested & (result["adj_p"] < 0.05) & (abs(result["haplotype1_frac"] - 0.5) >= 0.2), "status"
    ] = "significant_imbalance"
    return result[OUTPUT_COLUMNS]


def run(args: argparse.Namespace) -> int:
    ase = read_ase(args.ase_counts)
    phase = read_haplotype_map(args.haplotype_map)
    peaks = read_peaks(args.peaks)
    if ase.empty:
        empty_output(args.output)
        print("Evaluated peak/phase sets: 0 (ASE input contains no rows)")
        return 0
    ase = ase.loc[ase["totalCount"] >= 30].copy()
    if ase.empty:
        empty_output(args.output)
        print("Evaluated peak/phase sets: 0 (all sites are below totalCount 30)")
        return 0
    keys = ["contig", "position", "refAllele", "altAllele"]
    ase = ase.merge(
        phase[keys + ["phaseSet", "haplotype1Allele"]], on=keys,
        how="left", validate="one_to_one",
    )
    if ase[["phaseSet", "haplotype1Allele"]].isna().any().any():
        missing = ase.loc[ase["phaseSet"].isna(), keys]
        preview = ", ".join(
            f"{row.contig}:{row.position}:{row.refAllele}>{row.altAllele}"
            for row in missing.head(3).itertuples()
        )
        raise InputContractError(f"haplotype map is missing depth-eligible ASE variants: {preview}")
    ase["haplotype1Count"] = np.where(
        ase["haplotype1Allele"] == "REF", ase["refCount"], ase["altCount"]
    )
    ase["haplotype2Count"] = np.where(
        ase["haplotype1Allele"] == "REF", ase["altCount"], ase["refCount"]
    )
    mapped = map_variants_to_peaks(ase, peaks)
    if mapped.empty:
        empty_output(args.output)
        print("Evaluated peak/phase sets: 0 (no depth-eligible sites overlap peaks)")
        return 0
    result = evaluate_groups(mapped)
    result.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    print(
        f"Evaluated peak/phase sets: {len(result)}; "
        f"significant imbalance: {(result['status'] == 'significant_imbalance').sum()}"
    )
    return 0


def main() -> int:
    args = parse_args()
    try:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        return run(args)
    except (InputContractError, OSError, pd.errors.ParserError) as exc:
        print(f"input contract error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
