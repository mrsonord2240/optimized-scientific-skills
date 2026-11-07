# HTML and Cytoscape

```bash
python examples/interactive_network.py --graphml network.graphml --output-dir out/interactive
```

This creates basic and styled PyVis HTML and uses a copy for conversion to avoid mutating edge attributes. Direction is enabled when the input is directed, but the final community summary calls the undirected algorithm on the original graph and may fail for directed input; adapt or omit that summary and inspect arrowheads. Around 2,000 nodes is a practical starting point, not a guarantee of browser responsiveness. For desktop export, start Cytoscape and use `python examples/cytoscape_automation.py --graphml network.graphml --output-dir out/cytoscape`. It creates and styles a network in the active session and overwrites named exports; its PPI style suppresses arrows, so directed graphs require explicit arrow and sign mappings. See [details](../references/interactive-and-cytoscape.md).

Done when direction, graph attributes, styling, browser responsiveness or exported files are checked in the actual output.
