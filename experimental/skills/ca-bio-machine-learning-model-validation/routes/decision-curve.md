# Decision curve: net benefit

```bash
python scripts/decision_curve.py predictions.csv --label label --prob prob --thresholds 0.05,0.1,0.2,0.3,0.4,0.5
```

`predictions.csv` as in `routes/calibration.md`: held-out labels and predicted probabilities. The script prints net benefit of the model, treat-all and treat-none at each threshold. Net benefit = TP/n - (FP/n) x pt/(1-pt), where `pt` is the threshold probability encoding how many false positives one true positive is worth.

- The model is useful only at thresholds where its net benefit exceeds both treat-all and treat-none.
- Choose thresholds from clinical reasoning before looking at the output; do not pick the threshold that looks best.
- The curve assumes calibrated probabilities. Run `routes/calibration.md` first; recalibrate if the slope is far from 1.
- High AUC with no net benefit over treat-all at every plausible threshold means the model is not useful.

Done when you have reported the thresholds used, the net benefit of the model and of treat-all, and where the model beats both.
