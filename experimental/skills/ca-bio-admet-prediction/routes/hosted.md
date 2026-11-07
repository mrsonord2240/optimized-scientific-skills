# Hosted ADMETlab results

```python
from scripts.predict_admet import load_admetlab_results
results = load_admetlab_results('admetlab3_results.csv', task_id_column='task_id')
```

Obtain the current official API contract before submission and use its actual task-ID column. The helper validates a returned CSV; it does not submit a task. Modeled-endpoint uncertainty does not apply to calculated rules or properties. See `references/admetlab.md`.

Done when the submitted task and validated result are recorded.
