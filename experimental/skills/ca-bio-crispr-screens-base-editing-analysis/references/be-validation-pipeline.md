# Broad be-validation-pipeline

Read when post-processing CRISPResso2 BE amplicon output with the Broad notebooks.

## Broad be-validation-pipeline

The Broad Institute's `be-validation-pipeline` (https://broadinstitute.github.io/be-validation-pipeline/) is a CRISPResso2 post-processing and validation toolkit for BE amplicon data -- a set of Jupyter notebooks, not a workflow-engine pipeline. Run CRISPResso2 first, then execute the notebooks in order:

```bash
git clone https://github.com/broadinstitute/be-validation-pipeline
cd be-validation-pipeline
pip install -r requirements.txt

# Step 1: run CRISPResso2 in batch mode (or use the BEV tool on GPP LIMS).
# The batch file is tab-delimited with columns: name, fastq_r1, amplicon_seq, guide_seq
# (plus optional -w, -wc, --exclude_bp_from_left/right).
docker run -v ${PWD}:/DATA -w /DATA -i pinellolab/crispresso2 \
    CRISPRessoBatch --batch_settings batch_file.txt --skip_failed --base_editor_output

# Step 2: run the notebooks in order against the CRISPResso2 output
#   notebooks/01_BEV_allele_frequencies.ipynb
#   notebooks/02_BEV_nucleotide_percentage_plots.ipynb
#   notebooks/03_BEV_editing_efficiency.ipynb
# Outputs: allele-frequency tables, nucleotide-percentage plots, editing-efficiency heat maps
```

The notebooks cover allele-frequency tabulation, nucleotide-level editing quantification and editing-efficiency summaries. Hit calling is NOT part of this toolkit -- score the screen separately with drugZ or MAGeCK (drugZ is more sensitive than MAGeCK for drug-modifier chemogenomic screens like Hanna 2021 PARPi). Reuse the notebooks before writing custom BE amplicon parsers.

Required inputs for a screen: amplicon FASTQ (per sgRNA or pool), a library file (per-sgRNA spacer, target base, target amino acid, predicted bystander pattern), the BE chemistry (CBE or ABE), and ClinVar/COSMIC annotation for variant attribution.
