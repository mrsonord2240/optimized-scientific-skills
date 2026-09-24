# BE-Hive Editing-Efficiency Prediction

Read when predicting per-spacer editing efficiency or bystander outcomes with BE-Hive (Arbab 2020).

## BE-Hive Editing-Efficiency Prediction

BE-Hive (`be_predict_bystander`) takes a fixed **50nt substrate**, not the bare 20nt spacer:
19nt of upstream context + the 20nt spacer + the 3nt PAM + 8nt of downstream context
(`-19..30` in the README's own numbering; spacer occupies substrate positions 1-20 of that
window, PAM at 21-23). Getting the upstream-context length wrong by even 1nt silently
shifts every predicted position by one base with no error -- checked on this Skill's own
worked example, verified via BE-Hive's own `pred_df` diagnostic field.

```bash
python scripts/behive_predict.py --behive-parent /path/to/clone-parent-dir \
    --spacer TGATCACGTAGCATGCACGT --pam TGG --upstream ATGCATGGATCGTAGCTAG --downstream CATGCTAG \
    --editor BE4 --celltype mES   # run in the BE-Hive Python env
```

`scripts/behive_predict.py` builds the 50nt substrate, asserts the spacer offset (`substrate[19:39] == spacer`) before trusting `pred_df`'s position-labeled columns (e.g. `C4`, `C6`), and prints `Total predicted probability` and the top outcomes. A wrong substrate length or offset produces a plausible-looking but silently mis-positioned prediction.

Checked on BE-Hive git HEAD (maxwshen/be_predict_bystander, 2026-09-16) against a synthetic guide
with a known target C at spacer position 5 and bystander C at spacer position 7: `pred_df`'s `C4`/`C6`
columns (BE-Hive's own 0-indexed-from-position-4 editable-C naming) correctly identified both, and
`Total predicted probability` was 0.97-0.98 (not a stub, not all-zero).
