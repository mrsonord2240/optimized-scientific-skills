# Installing the tools

```bash
conda install -c bioconda mageck mageck-vispr      # not on PyPI
git clone https://github.com/hart-lab/bagel        # BAGEL2: BAGEL.py, CEGv2.txt, NEGv1.txt
git clone https://github.com/hart-lab/drugz
git clone https://github.com/felicityallen/JACKS   # the PyPI `jacks` is unrelated
pip install crispr-chronos                         # import chronos; TensorFlow, own environment
pip install pandas numpy matplotlib scipy
R -e "devtools::install_github('francescojm/CRISPRcleanR')"
R -e "remotes::install_github('WubingZhang/MAGeCKFlute')"   # removed from Bioconductor at 3.22
```

Check versions with `mageck --version`, `BAGEL.py fc --help`, `drugz.py -h`, `packageVersion('CRISPRcleanR')`. If a call raises ImportError, AttributeError or TypeError, introspect the installed package and adapt rather than retrying.
