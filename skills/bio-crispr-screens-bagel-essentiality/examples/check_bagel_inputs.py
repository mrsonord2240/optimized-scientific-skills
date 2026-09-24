#!/usr/bin/env python3
"""Pre-flight and post-run sanity checks for BAGEL.py bf.

BAGEL.py gives no runtime signal for the two silent/ugly reference-set mistakes:
  - species / symbol mismatch  -> raw scipy traceback (gaussian_kde) or too few reference genes
  - -e / -n swapped            -> exit 0 and a bayes_factor.txt whose BF column is all `nan`

Usage (checked on BAGEL2 build 115, pandas 2.x/3.x):
  check_bagel_inputs.py pre  foldchange.foldchange CEGv2.txt NEGv1.txt Sample1,Sample2,Sample3
  check_bagel_inputs.py post bayes_factor.txt

Exit code 1 = do not trust / do not run; 0 = passed (warnings may still print).
"""
import sys
import pandas as pd

MIN_REF_OVERLAP = 100      # BAGEL needs enough reference genes present in the screen for a KDE
MAX_NAN_FRAC = 0.05        # post-run: more than this fraction of NaN BF is a failed run


def read_genes(path):
    """Reference files are tab-delimited with a GENE header (CEGv2.txt has CRLF endings)."""
    df = pd.read_csv(path, sep='\t', dtype=str)
    return set(df.iloc[:, 0].str.strip())


def preflight(fc_path, ceg_path, neg_path, treatment):
    fc = pd.read_csv(fc_path, sep='\t')
    gene_col = 'GENE'
    missing = [c for c in treatment if c not in fc.columns]
    if missing:
        print(f"FAIL: treatment column(s) {missing} not in {fc_path}; columns are {list(fc.columns[2:])}")
        return 1
    ceg, neg = read_genes(ceg_path), read_genes(neg_path)
    screen_genes = set(fc[gene_col].astype(str))
    n_ceg, n_neg = len(ceg & screen_genes), len(neg & screen_genes)
    print(f"reference overlap with screen GENE column: essential {n_ceg}/{len(ceg)}, non-essential {n_neg}/{len(neg)}")
    bad = 0
    for label, n, tot in (('essential (-e)', n_ceg, len(ceg)), ('non-essential (-n)', n_neg, len(neg))):
        if n < MIN_REF_OVERLAP:
            print(f"FAIL: only {n}/{tot} {label} genes are in the screen -- species or gene-symbol mismatch "
                  "(BAGEL.py would crash with 'dataset input should have multiple elements' or give a useless KDE)")
            bad = 1
    if bad:
        return 1
    # Dropout sanity: essentials must have a lower mean LFC than non-essentials.
    lfc = fc[treatment].mean(axis=1)
    m_ess = lfc[fc[gene_col].isin(ceg)].mean()
    m_non = lfc[fc[gene_col].isin(neg)].mean()
    print(f"mean LFC: -e genes {m_ess:.3f}, -n genes {m_non:.3f}")
    if m_ess >= m_non:
        print("FAIL: -e genes do not drop out more than -n genes. Either -e/-n are swapped (BAGEL.py would "
              "write an all-`nan` BF column and exit 0) or the screen shows no essential-gene dropout "
              "(wrong -c columns, control/treatment reversed, or a failed screen).")
        return 1
    print("pre-flight OK")
    return 0


def postrun(bf_path):
    df = pd.read_csv(bf_path, sep='\t')
    # BAGEL.py writes nan with a leading space (' nan'), so the column can load as text
    bf = pd.to_numeric(df['BF'], errors='coerce')
    frac = bf.isna().mean()
    print(f"{bf.notna().sum()}/{len(bf)} genes have a numeric BF ({frac:.1%} NaN)")
    if frac > MAX_NAN_FRAC:
        print("FAIL: BF column is mostly NaN. Likely causes: -e/-n swapped, species/gene-symbol mismatch "
              "between the reference files and the screen, or a wrong -c treatment column. "
              "Do not read this as 'no essential genes'.")
        return 1
    print(f"BF>6: {(bf > 6).sum()} genes; median BF {bf.median():.2f}")
    if bf.median() > -1 and (bf > 6).sum() == 0:
        print("WARN: no gene has BF>6 -- thin/unrepresentative reference set, or a genuinely flat screen.")
    print("post-run OK")
    return 0


if __name__ == '__main__':
    if len(sys.argv) == 6 and sys.argv[1] == 'pre':
        sys.exit(preflight(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5].split(',')))
    if len(sys.argv) == 3 and sys.argv[1] == 'post':
        sys.exit(postrun(sys.argv[2]))
    print(__doc__)
    sys.exit(2)
