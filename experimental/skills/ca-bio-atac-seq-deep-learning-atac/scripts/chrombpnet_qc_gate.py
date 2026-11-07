#!/usr/bin/env python3
"""QC gate for a trained chromBPNet model. Exit 0 = pass, 1 = fail, 2 = QC file missing.

Thresholds are the ones chromBPNet 1.0.1 prints in its own QC report
(chrombpnet/helpers/generate_reports/make_html.py, marginal_footprinting.py):
  bias model, pearsonr in peaks   > -0.3   (chrombpnet aborts below -0.5)
  corrected model, pearsonr peaks > 0.5
  Tn5 marginal footprint, max     < 0.003  (file starts with 'corrected_')

usage: chrombpnet_qc_gate.py BIAS_METRICS.json CHROMBPNET_METRICS.json MAX_BIAS_RESPONSE.txt
"""
import json
import sys


def main(bias_json, model_json, response_txt):
    try:
        bias_r = json.load(open(bias_json, encoding="utf-8"))["counts_metrics"]["peaks"]["pearsonr"]
        model = json.load(open(model_json, encoding="utf-8"))
        # train writes the key "peaks", pred_bw -bw writes "regions"
        counts, prof = model["counts_metrics"], model["profile_metrics"]
        key = "peaks" if "peaks" in counts else "regions"
        model_r = counts[key]["pearsonr"]
        jsd = prof[key]["median_norm_jsd"]
        response = open(response_txt, encoding="utf-8").read().strip()
    except (OSError, KeyError, ValueError) as exc:
        print("QC files missing or unreadable: %r" % exc)
        return 2
    checks = [
        ("bias model pearsonr in peaks > -0.3", bias_r, bias_r > -0.3),
        ("corrected model pearsonr in peaks > 0.5", model_r, model_r > 0.5),
        ("Tn5 footprint max response < 0.003", response, response.startswith("corrected_")),
    ]
    for name, value, ok in checks:
        print("%s  %s  (%s)" % ("PASS" if ok else "FAIL", name, value))
    print("INFO  corrected model median_norm_jsd in peaks = %.3f (higher is better; depth sensitive)" % jsd)
    return 0 if all(ok for _, _, ok in checks) else 1


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
