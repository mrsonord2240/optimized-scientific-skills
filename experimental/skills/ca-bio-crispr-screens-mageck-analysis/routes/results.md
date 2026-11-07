# Summarize MLE and plot RRA
~~~bash
python scripts/time_course.py timecourse_mle.gene_summary.txt time_course.tsv --conditions day7,day14,day21
~~~
The source function selects monotonic beta trends with any FDR <0.05; despite computing same-sign flags it does not require them. Review sign and number of conditions separately. For an existing RRA summary:
~~~bash
python scripts/volcano.py drug_vs_veh.gene_summary.txt volcano.png --direction neg
~~~
The source plotting defaults are FDR 0.05 and absolute LFC 1. MAGeCKFlute and VISPR integration examples and limitations are in [visualization notes](../references/visualization.md); verify package versions and supply context-appropriate parameters.
Done when summaries, FDR behavior, plotting assumptions and limitations are reviewed.
