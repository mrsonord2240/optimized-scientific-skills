"""Bounded ACMG/AMP evidence helpers for training and method review.

These functions are not a diagnostic classifier. They expose selected published
thresholds for an evidence ledger that still requires qualified review.
"""
from __future__ import annotations

import math
import re
from collections.abc import Sequence
from datetime import date
from numbers import Real
from typing import Any

import requests

NON_DIAGNOSTIC_NOTICE = (
    "Educational evidence sketch only; not for patient diagnosis, treatment, or "
    "stand-alone clinical reporting. A qualified variant scientist or clinical "
    "geneticist must review disease, transcript, VCEP rules, provenance, conflicts, and uncertainty."
)
PVS1_REVIEW_REQUIRED = "PVS1_REVIEW_REQUIRED"

REVEL_THRESHOLDS = [
    (float("-inf"), 0.003, True, True, "BP4_VeryStrong"),
    (0.003, 0.016, False, True, "BP4_Strong"),
    (0.016, 0.183, False, True, "BP4_Moderate"),
    (0.183, 0.290, False, True, "BP4_Supporting"),
    (0.290, 0.644, False, False, None),
    (0.644, 0.773, True, False, "PP3_Supporting"),
    (0.773, 0.932, True, False, "PP3_Moderate"),
    (0.932, float("inf"), True, True, "PP3_Strong"),
]
BAYESDEL_NOAF_THRESHOLDS = [
    (float("-inf"), -0.36, True, True, "BP4_Moderate"),
    (-0.36, -0.18, False, True, "BP4_Supporting"),
    (-0.18, 0.13, False, False, None),
    (0.13, 0.27, True, False, "PP3_Supporting"),
    (0.27, 0.50, True, False, "PP3_Moderate"),
    (0.50, float("inf"), True, True, "PP3_Strong"),
]

# Bergquist's +/-3 intervals are explicit rather than rounded to a legacy label.
STRENGTH_POINTS = {
    "PVS1_VeryStrong": 8, "PVS1_Strong": 4, "PVS1_Moderate": 2, "PVS1_Supporting": 1,
    "PS1": 4, "PS2": 4, "PS3_VeryStrong": 8, "PS3": 4,
    "PS3_Moderate": 2, "PS3_Supporting": 1, "PS4": 4,
    "PS4_Moderate": 2, "PS4_Supporting": 1, "PM1": 2,
    "PM2_Supporting": 1, "PM2": 1, "PM3": 2, "PM3_Strong": 4,
    "PM3_VeryStrong": 8, "PM4": 2, "PM5": 2, "PM6": 2,
    "PP1": 1, "PP1_Moderate": 2, "PP1_Strong": 4, "PP2": 1,
    "PP3_Supporting": 1, "PP3_Moderate": 2, "PP3_3pt": 3,
    "PP3_Strong": 4, "PP4": 1, "BA1": -100,
    "BS1": -4, "BS2": -4, "BS3": -4, "BS3_Moderate": -2,
    "BS3_Supporting": -1, "BS4": -4, "BP1": -1, "BP2": -1,
    "BP3": -1, "BP4_Supporting": -1, "BP4_Moderate": -2,
    "BP4_3pt": -3, "BP4_Strong": -4, "BP4_VeryStrong": -8,
    "BP5": -1, "BP7": -1,
}
RETIRED_CRITERIA = {"PP5", "BP6"}
EVIDENCE_FAMILIES = {code: code.split("_", 1)[0] for code in STRENGTH_POINTS}


def _real(name: str, value: Any) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def _probability(name: str, value: Any) -> float:
    value = _real(name, value)
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be between 0 and 1 inclusive")
    return value


def _boolean(name: str, value: Any) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be bool")
    return value


def pejaver_calibrate(score, thresholds):
    """Map a finite score through an explicit calibrated interval table."""
    if score is None:
        return None
    score = _real("score", score)
    for interval in thresholds:
        if len(interval) == 3:
            lo, hi, code = interval
            include_lo, include_hi = True, False
        elif len(interval) == 5:
            lo, hi, include_lo, include_hi, code = interval
        else:
            raise ValueError("threshold intervals must have 3 or 5 fields")
        above_lo = score > lo or (include_lo and score == lo)
        below_hi = score < hi or (include_hi and score == hi)
        if above_lo and below_hi:
            return code
    return None


def revel_pp3_bp4(revel_score):
    if revel_score is None:
        return None
    return pejaver_calibrate(_probability("revel_score", revel_score), REVEL_THRESHOLDS)


def bayesdel_pp3_bp4(bayesdel_noaf_score):
    return pejaver_calibrate(bayesdel_noaf_score, BAYESDEL_NOAF_THRESHOLDS)


def spliceai_walker2023(ds_max):
    if ds_max is None:
        return None
    ds_max = _probability("ds_max", ds_max)
    if ds_max >= 0.20:
        return "PP3_Supporting"
    if ds_max <= 0.10:
        return "BP4_Supporting"
    return None


def alphamissense_pp3_bp4(am_score):
    """Bergquist 2025 AlphaMissense calibration, including +/-3-point bands."""
    if am_score is None:
        return None
    score = _probability("am_score", am_score)
    if score <= 0.070:
        return "BP4_3pt"
    if score <= 0.099:
        return "BP4_Moderate"
    if score <= 0.169:
        return "BP4_Supporting"
    if score <= 0.791:
        return None
    if score <= 0.905:
        return "PP3_Supporting"
    if score <= 0.971:
        return "PP3_Moderate"
    if score <= 0.989:
        return "PP3_3pt"
    return "PP3_Strong"


# Compatibility name; behavior is now the complete calibrated mapping.
alphamissense_supporting_only = alphamissense_pp3_bp4


def pvs1_decision_tree(
    variant_type, is_nmd_predicted, coding_pct_removed, in_critical_region,
    is_disease_relevant_transcript, lof_mechanism_established, *,
    splice_consequence=None, rescue_transcript_excluded=None,
):
    """PVS1 sketch with mandatory mechanism, transcript, and splice-review gates.

    ``PVS1_REVIEW_REQUIRED`` is a stop state, not an evidence code.
    """
    allowed = {"nonsense", "frameshift", "splice_donor", "splice_acceptor",
               "initiation", "single_exon_del", "multi_exon_del"}
    if not isinstance(variant_type, str) or variant_type not in allowed:
        raise ValueError(f"variant_type must be one of {sorted(allowed)}")
    _boolean("is_nmd_predicted", is_nmd_predicted)
    removed = _probability("coding_pct_removed", coding_pct_removed)
    _boolean("in_critical_region", in_critical_region)
    _boolean("is_disease_relevant_transcript", is_disease_relevant_transcript)
    _boolean("lof_mechanism_established", lof_mechanism_established)
    if not lof_mechanism_established or not is_disease_relevant_transcript:
        return None

    if variant_type in {"splice_donor", "splice_acceptor"}:
        if splice_consequence not in {"nmd", "in_frame_disruptive", "no_impact"}:
            return PVS1_REVIEW_REQUIRED
        if is_nmd_predicted != (splice_consequence == "nmd"):
            raise ValueError(
                "is_nmd_predicted conflicts with splice_consequence; reconcile NMD review before PVS1"
            )
        if rescue_transcript_excluded is not True:
            if rescue_transcript_excluded is not None:
                _boolean("rescue_transcript_excluded", rescue_transcript_excluded)
            return PVS1_REVIEW_REQUIRED
        if splice_consequence == "no_impact":
            return None
        if splice_consequence == "nmd":
            return "PVS1_VeryStrong"
        return "PVS1_Strong" if in_critical_region or removed > 0.10 else "PVS1_Moderate"

    if splice_consequence is not None or rescue_transcript_excluded is not None:
        raise ValueError("splice review fields are only valid for canonical splice variants")
    if variant_type in {"nonsense", "frameshift"}:
        if is_nmd_predicted:
            return "PVS1_VeryStrong"
        return "PVS1_Strong" if in_critical_region or removed > 0.10 else "PVS1_Moderate"
    if variant_type == "initiation":
        return "PVS1_Moderate"
    return "PVS1_VeryStrong" if in_critical_region else "PVS1_Strong"


def whiffin_max_credible_af(prevalence, max_allelic_contribution=1.0,
                            max_genetic_contribution=1.0, penetrance=1.0):
    prevalence = _probability("prevalence", prevalence)
    allelic = _probability("max_allelic_contribution", max_allelic_contribution)
    genetic = _probability("max_genetic_contribution", max_genetic_contribution)
    penetrance = _probability("penetrance", penetrance)
    if penetrance == 0:
        raise ValueError("penetrance must be greater than zero")
    return (prevalence * genetic * allelic) / (penetrance * 2)


def bs1_ba1(grpmax_faf95, max_credible_af, ba1_threshold=0.05):
    maximum = _probability("max_credible_af", max_credible_af)
    ba1 = _probability("ba1_threshold", ba1_threshold)
    if ba1 == 0 or maximum > ba1:
        raise ValueError("require 0 < ba1_threshold and max_credible_af <= ba1_threshold")
    if grpmax_faf95 is None:
        return "PM2_Supporting"
    observed = _probability("grpmax_faf95", grpmax_faf95)
    if observed == 0:
        return "PM2_Supporting"
    if observed > ba1:
        return "BA1"
    if observed > maximum:
        return "BS1"
    return None


def ps3_oddspath(odds_path):
    """Exact Brnich 2020 Table 3 boundaries."""
    if odds_path is None:
        return None
    value = _real("odds_path", odds_path)
    if value <= 0:
        raise ValueError("odds_path must be greater than zero")
    if value > 350:
        return "PS3_VeryStrong"
    if value > 18.7:
        return "PS3"
    if value > 4.3:
        return "PS3_Moderate"
    if value > 2.1:
        return "PS3_Supporting"
    if value < 0.053:
        return "BS3"
    if value < 0.23:
        return "BS3_Moderate"
    if value < 0.48:
        return "BS3_Supporting"
    return None


def _criteria(criteria_assigned: Sequence[str]) -> list[str]:
    if isinstance(criteria_assigned, (str, bytes)) or not isinstance(criteria_assigned, Sequence):
        raise TypeError("criteria_assigned must be a sequence of evidence-code strings")
    criteria = list(criteria_assigned)
    if any(not isinstance(code, str) for code in criteria):
        raise TypeError("every evidence code must be a string")
    retired = sorted(set(criteria) & RETIRED_CRITERIA)
    if retired:
        raise ValueError(f"retired current-rule evidence code(s): {', '.join(retired)}")
    unknown = sorted(set(criteria) - STRENGTH_POINTS.keys())
    if unknown:
        raise ValueError(f"unknown evidence code(s): {', '.join(unknown)}")
    duplicates = sorted({code for code in criteria if criteria.count(code) > 1})
    if duplicates:
        raise ValueError(f"duplicate evidence code(s): {', '.join(duplicates)}")
    families: dict[str, list[str]] = {}
    for code in criteria:
        families.setdefault(EVIDENCE_FAMILIES[code], []).append(code)
    repeated_families = {family: codes for family, codes in families.items() if len(codes) > 1}
    if repeated_families:
        details = "; ".join(
            f"{family}: {', '.join(codes)}" for family, codes in sorted(repeated_families.items())
        )
        raise ValueError(f"multiple strengths or aliases from one evidence family: {details}")
    if "PP3" in families and "BP4" in families:
        raise ValueError("opposing PP3 and BP4 computational evidence requires conflict review")
    return criteria


def tavtigian_classify(criteria_assigned):
    """Sum validated points and return a result carrying a clinical-use guard."""
    criteria = _criteria(criteria_assigned)
    points = sum(STRENGTH_POINTS[code] for code in criteria)
    if "BA1" in criteria:
        category = "Benign"
    elif points >= 10:
        category = "Pathogenic"
    elif points >= 6:
        category = "Likely Pathogenic"
    elif points >= 0:
        category = "VUS"
    elif points >= -6:
        category = "Likely Benign"
    else:
        category = "Benign"
    return {"classification": category, "points": points, "criteria": criteria,
            "review_required": True, "clinical_use": NON_DIAGNOSTIC_NOTICE}


def classify_with_subsumption(criteria_assigned):
    criteria = _criteria(criteria_assigned)
    if any(code.startswith("PVS1") for code in criteria):
        criteria = [c for c in criteria if not c.startswith("PP3") and c != "PM4"]
    return tavtigian_classify(criteria)


def genebe_api(chrom, pos, ref, alt, genome="hg38", timeout=30):
    """Current GeneBe public single-variant coordinate interface."""
    if not isinstance(chrom, str) or not re.fullmatch(r"(?:chr)?(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)", chrom):
        raise ValueError("chrom must be a canonical chromosome name")
    if isinstance(pos, bool) or not isinstance(pos, int) or pos <= 0:
        raise ValueError("pos must be a positive integer")
    for name, allele in (("ref", ref), ("alt", alt)):
        if not isinstance(allele, str) or not re.fullmatch(r"[ACGTN]+", allele.upper()):
            raise ValueError(f"{name} must contain one or more DNA bases")
    if genome not in {"hg19", "hg38"}:
        raise ValueError("genome must be 'hg19' or 'hg38'")
    timeout = _real("timeout", timeout)
    if timeout <= 0:
        raise ValueError("timeout must be greater than zero")
    response = requests.get("https://api.genebe.net/cloud/api-public/v1/variant",
                            params={"chr": chrom.removeprefix("chr"), "pos": pos,
                                    "ref": ref.upper(), "alt": alt.upper(), "genome": genome},
                            timeout=timeout)
    response.raise_for_status()
    return response.json()


def cspec_gene_versions(gene_symbol, timeout=30):
    """Fetch current version records for a ClinGen CSpec gene identifier."""
    if not isinstance(gene_symbol, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.-]*", gene_symbol):
        raise ValueError("gene_symbol must be a nonempty HGNC-style symbol")
    timeout = _real("timeout", timeout)
    if timeout <= 0:
        raise ValueError("timeout must be greater than zero")
    url = ("https://cspec.genome.network/cspec/Gene/id/" + gene_symbol.upper()
           + "/SequenceVariantInterpretation/version")
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.json()


def cancer_amp_tier(evidence_level, *, same_tumor_type, clinical_significance,
                    knowledgebase, access_date):
    """Assign a bounded AMP tier from explicit evidence and current provenance."""
    levels = {"regulatory_approved", "professional_guideline", "clinical_trial", "preclinical", "none"}
    significance = {"oncogenic", "likely_oncogenic", "uncertain", "likely_benign", "benign"}
    if evidence_level not in levels:
        raise ValueError(f"evidence_level must be one of {sorted(levels)}")
    _boolean("same_tumor_type", same_tumor_type)
    if clinical_significance not in significance:
        raise ValueError(f"clinical_significance must be one of {sorted(significance)}")
    if not isinstance(knowledgebase, str) or not knowledgebase.strip():
        raise ValueError("knowledgebase and its evidence record are required")
    if not isinstance(access_date, str) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", access_date):
        raise ValueError("access_date must use YYYY-MM-DD")
    try:
        parsed_access_date = date.fromisoformat(access_date)
    except ValueError as exc:
        raise ValueError("access_date must be a real ISO calendar date") from exc
    if parsed_access_date.isoformat() != access_date:
        raise ValueError("access_date must use canonical YYYY-MM-DD")
    if evidence_level != "none" and clinical_significance in {"likely_benign", "benign"}:
        raise ValueError("actionable evidence conflicts with benign significance; review required")
    if evidence_level == "regulatory_approved" and same_tumor_type:
        tier, rationale = "Tier I-A", "regulatory approval for this biomarker and tumor type"
    elif evidence_level == "professional_guideline" and same_tumor_type:
        tier, rationale = "Tier I-B", "professional guideline evidence for this tumor type"
    elif evidence_level in {"regulatory_approved", "professional_guideline", "clinical_trial"}:
        tier, rationale = "Tier II-C", "different-tumor or investigational clinical evidence"
    elif evidence_level == "preclinical":
        tier, rationale = "Tier II-D", "preclinical evidence only"
    elif clinical_significance == "uncertain":
        tier, rationale = "Tier III", "uncertain clinical significance"
    elif clinical_significance in {"likely_benign", "benign"}:
        tier, rationale = "Tier IV", "benign or likely benign"
    else:
        raise ValueError("oncogenic significance without evidence requires review, not a tier")
    return {"tier": tier, "rationale": rationale, "knowledgebase": knowledgebase.strip(),
            "access_date": access_date, "review_required": True,
            "clinical_use": NON_DIAGNOSTIC_NOTICE}


if __name__ == "__main__":
    # One computational predictor only; context intentionally unresolved.
    criteria = [revel_pp3_bp4(0.95), ps3_oddspath(8.0), bs1_ba1(0.0, 1e-6)]
    retained = [code for code in criteria if code is not None]
    result = classify_with_subsumption(retained)
    print("NON-DIAGNOSTIC TRAINING EXAMPLE")
    print(NON_DIAGNOSTIC_NOTICE)
    print("Context: disease=unresolved; transcript=unresolved; VCEP=not checked; assay validity=unreviewed")
    print(f"Criteria: {retained}")
    print(f"Illustrative points/category: {result['points']} / {result['classification']}")
    print("Uncertainty: do not use until context and conflicts receive qualified review.")
