# Reference: drugZ Aug-2019+ (hart-lab/drugz), pandas 2.2+, numpy 1.26+ | Verify API if version differs
#
# drugZ chemogenomic screen analysis end-to-end.
# Compares drug-treated vs vehicle-treated arms at matched timepoint.

import subprocess
import sys

import pandas as pd
import numpy as np
from pathlib import Path

# === INPUTS ===
counts_file = 'counts.txt'                          # sgRNA, GENE, sample columns
vehicle_samples = 'Veh_r1,Veh_r2,Veh_r3'
drug_samples = 'Drug_r1,Drug_r2,Drug_r3'
output_dir = Path('drugz_output')
output_dir.mkdir(exist_ok=True)

# === DOWNLOAD CEGv2 (for excluding essentials from null) ===
ceg_file = Path('CEGv2.txt')
if not ceg_file.exists():
    subprocess.run(['curl', '-L', '-o', str(ceg_file),
                    'https://raw.githubusercontent.com/hart-lab/bagel/master/CEGv2.txt'],
                   check=True)

# === STEP 1: STANDARD drugZ RUN ===
output_file = output_dir / 'drugz_standard.txt'
subprocess.run([sys.executable, 'drugz.py',   # the interpreter running this script, not whatever 'python' resolves to
                '-i', counts_file,
                '-o', str(output_file),
                '-c', vehicle_samples,
                '-x', drug_samples,
                '-p', '5'], check=True)

# === STEP 2: drugZ EXCLUDING ESSENTIALS FROM THE ANALYSIS ===
# drugZ's -r takes a COMMA-DELIMITED GENE LIST, not a file path, and it drops those genes'
# guides before Z-scoring, so they also vanish from the output. CEGv2.txt is tab-separated
# with a header (GENE, HGNC_ID, ENTREZ_ID): take column 1, skip the header. Joining the raw
# lines produces tokens like "AARS\tHGNC:20\t16", which match nothing and silently exclude
# nothing while the tool still exits 0.
ceg_genes_list = [line.split('\t')[0].strip()
                  for line in ceg_file.read_text().splitlines()[1:] if line.strip()]
ceg_genes = ','.join(ceg_genes_list)
output_file_clean = output_dir / 'drugz_ceg_excluded.txt'
subprocess.run([sys.executable, 'drugz.py',   # the interpreter running this script, not whatever 'python' resolves to
                '-i', counts_file,
                '-o', str(output_file_clean),
                '-c', vehicle_samples,
                '-x', drug_samples,
                '-r', ceg_genes,
                '-p', '5'], check=True)

# Confirm the exclusion actually happened -- drugZ reports nothing when -r matches no gene.
genes_all = set(pd.read_csv(output_file, sep='\t')['GENE'])
genes_clean = set(pd.read_csv(output_file_clean, sep='\t')['GENE'])
removed = genes_all - genes_clean
print(f'CEGv2 genes excluded from the analysis: {len(removed)} of {len(ceg_genes_list)} listed')
assert removed, '-r matched no gene: check the reference-file parsing'

# === STEP 3: INTERPRET ===
df = pd.read_csv(output_file, sep='\t')

# fdr_synth = sensitizer (negative Z; gene KO sensitizes to drug)
# fdr_supp = suppressor (positive Z; gene KO confers resistance)
sensitizers = df[df['fdr_synth'] < 0.05].sort_values('normZ').head(50)
suppressors = df[df['fdr_supp'] < 0.05].sort_values('normZ', ascending=False).head(50)

# Save tier-stratified output
sensitizers.to_csv(output_dir / 'top50_sensitizers.tsv', sep='\t', index=False)
suppressors.to_csv(output_dir / 'top50_suppressors.tsv', sep='\t', index=False)

# === STEP 4: SUMMARY ===
print(f'Sensitizers (drug + KO synergistic):  {len(sensitizers)}')
print(f'Suppressors (KO confers resistance):  {len(suppressors)}')
print(f'Top 5 sensitizers: {", ".join(sensitizers.head(5)["GENE"].astype(str))}')
print(f'Top 5 suppressors: {", ".join(suppressors.head(5)["GENE"].astype(str))}')

# === STEP 5: HIGH-CONFIDENCE HITS ===
# Apply stringent threshold: fdr_synth < 0.01 AND normZ < -3
strong_sens = df[(df['fdr_synth'] < 0.01) & (df['normZ'] < -3)]
print(f'Strong sensitizers (fdr<0.01, normZ<-3): {len(strong_sens)}')

# === STEP 6: DOSE CONSISTENCY (if multi-dose) ===
# Run drugZ once per dose against the same vehicle, then keep genes with one sign across
# doses that also pass FDR at the top dose. See SKILL.md's dose-consistency rule.

def per_dose_drugz(doses_dict, top_dose, fdr=0.05, direction='synth'):
    """doses_dict: {'low': 'Drug_low_r1,Drug_low_r2', 'mid': ..., 'high': ...}"""
    dose_files = {}
    for dose, samples in doses_dict.items():
        out = output_dir / f'drugz_{dose}.txt'
        subprocess.run([sys.executable, 'drugz.py',   # the interpreter running this script, not whatever 'python' resolves to
                        '-i', counts_file,
                        '-o', str(out),
                        '-c', vehicle_samples,
                        '-x', samples,
                        '-p', '5'], check=True)
        dose_files[dose] = out

    frames = {d: pd.read_csv(f, sep='\t').set_index('GENE') for d, f in dose_files.items()}
    normz = pd.DataFrame({d: f['normZ'] for d, f in frames.items()}).dropna()
    same_sign = normz.gt(0).all(axis=1) | normz.lt(0).all(axis=1)
    passes_top = (frames[top_dose][f'fdr_{direction}'] < fdr).reindex(normz.index).fillna(False)
    hits = normz[same_sign & passes_top].copy()
    hits['normZ_top_dose'] = frames[top_dose]['normZ'].reindex(hits.index)
    return hits.sort_values('normZ_top_dose')

# hits = per_dose_drugz({'low': 'Drug_low_r1,Drug_low_r2',
#                        'mid': 'Drug_mid_r1,Drug_mid_r2',
#                        'high': 'Drug_high_r1,Drug_high_r2'}, top_dose='high')
