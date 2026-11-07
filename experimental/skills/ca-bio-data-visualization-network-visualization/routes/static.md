# Static network layouts

```bash
python examples/network_plots.py --graphml network.graphml --output-dir out/static
```

Use an attributed undirected GraphML graph. The source draws degree, communities, hubs and score views with seed 42; missing weights default to 1 and widths scale to 0.5–4 points. Score thresholds are 700 and 900; absent scores default to zero in the confidence view. Labels are degree-ranked with a cap of 30; add repeatable `--label NODE` for requested nodes. Existing output names may be overwritten. For comparisons run `python scripts/layouts.py network.graphml --output-dir out/layouts --layouts spring,kamada-kawai,circular,spectral,forceatlas2`; native ForceAtlas2 requires NetworkX 3.4+. For directed signed graphs use `python scripts/directed_network.py network.graphml --output out/directed.png --sign-attribute sign`, requiring Graphviz and pydot. Unknown explicit signs are rejected, but missing signs default to `+`; validate completeness before drawing. See [layout details](../references/python-layouts.md).

Done when graph size, direction, attributes, layout seed, scales and output are inspected and documented.
