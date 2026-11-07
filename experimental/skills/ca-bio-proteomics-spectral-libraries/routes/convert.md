# Convert library tables
~~~bash
python scripts/spectronaut_to_diann.py --in spectronaut.tsv --out diann.tsv
~~~
The source maps named columns and checks iRT lies between -50 and 200; a numerical range alone does not prove units or format compatibility. Verify modification notation, relative intensities, required columns and target engine version. For a source OpenSWATH transition template:
~~~bash
python scripts/build_openswath_tsv.py --peptides peptides.tsv --out library.tsv --n-frag 6
~~~
Inputs are sequence, charge, protein and irt. This computes theoretical y-ion m/z and assigns illustrative 1000/i intensities; it is not a measured or predicted intensity library. Conversion to TraML and decoy generation require additional validated tools and checks in [format conversion](../references/format-conversion.md).
Done when units, transition identities, intensities and decoy compatibility are reviewed.
