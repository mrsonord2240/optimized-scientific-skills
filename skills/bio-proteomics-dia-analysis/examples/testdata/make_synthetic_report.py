"""Write a tiny DIA-NN-shaped Parquet fixture for the SKILL.md filter example.

The fixture contains a group that passes the run-level q-values but fails only
the global protein-group q-value.  It is deliberately synthetic and is not a
mass-spectrometry result.
"""
from pathlib import Path

import pandas as pd


OUT = Path(__file__).with_name("synthetic_report.parquet")
rows = [
    {
        "Protein.Group": group,
        "Run": run,
        "Q.Value": 0.005,
        "PG.Q.Value": 0.006,
        "Global.Q.Value": 0.007,
        "Global.PG.Q.Value": global_pg_q,
        "PG.MaxLFQ": value,
    }
    for group, global_pg_q, values in [
        ("P001", 0.008, (1000.0, 1200.0)),
        ("P002", 0.009, (0.0, 800.0)),
        ("LOWCONF_GLOBAL_ONLY", 0.02, (900.0, 900.0)),
    ]
    for run, value in zip(("control", "treated"), values)
]

pd.DataFrame(rows).to_parquet(OUT, index=False)
print(OUT)
