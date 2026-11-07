# Several samples per patient or donor

```bash
python scripts/cv_auc.py data.csv --label label --group patient_id --drop sample_id
```

`--group` names the patient, donor, tumour or replicate column. Whole groups move between train and test together.

- When the label belongs to the patient (constant within a group), the script stratifies the patients, not the rows, and scores one AUC per patient by averaging that patient's scores. The output line says `scored per group`.
- When the label varies within a patient, it falls back to best-effort `StratifiedGroupKFold`; with few groups a test fold can be single-class.
- Do not leave `--group` out and rely on a random split. Rows from one patient in both train and test inflate the AUC.
- Tuning anything too: use `routes/nested-cv.md` with `--group`.

Done when you have reported mean AUC and SD from the script's output line, with the number of groups.
