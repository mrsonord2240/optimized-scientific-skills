# Provenance

This normalized Skill derives from:

- Repository: `GPTomics/bioSkills`
- Commit: `d91ed3d563019e649dc854c56ccd62551359488a`
- Source path: `clip-seq/ago-clip-mirna-targets`
- Repository license: MIT
- Copyright notice: `Copyright (c) 2026 GPTomics`
- Provider author metadata: GPTomics
- Last source-path commit author at the pinned revision: Domen Jemec

The optimized-skill repository carries the provider's MIT license and
copyright notice in its root `LICENSE`. Normalization preserves provider
attribution in `SKILL.md` frontmatter.

Canonical Git blob identities at the pinned revision are the portable source
identity. Checkout SHA-256 values can vary with newline conversion and are not
used as source Git identities.

| Source file | Git blob |
|---|---|
| `SKILL.md` | `0b35c982eb5aa65c575eb56b30e697073cb0423e` |
| `usage-guide.md` | `18d68dcf8b7b4b67318f53e87a282e23ae2c2f2e` |
| `examples/run_chimeric_eclip.sh` | `a28f9d5955b47aafa03a94da72113880d5a7c488` |

Normalization consolidates repeated method-selection and usage guidance into
routed references. The standalone shell example is retained as a local script;
its original blob identity remains recorded above. The optimization repair
rewrites that script to the current pinned Hyb interface and adds independently
testable parser and coordinate-conversion surfaces; it does not attribute those
repairs to the provider.

## Optimization interface evidence

- Hyb executable: `gkudla/hyb` public commit
  `028ab6371ce793ca5e86f475fce1f2cc6ad3c677` (GPL-3.0 source, invoked as an
  external tool and not redistributed here).
- Yeo chimeric-eCLIP interface classification: public MIT repository commit
  `75fe74e90e6e4ca670a5af76836d80db09bdbcb1`; total and targeted read-layout
  routes are kept distinct.
- HEAP: Li et al. 2020, DOI `10.1016/j.molcel.2020.05.009`, GEO `GSE139349`;
  associated CLIPanalyze source commit
  `fcab2db41b0c29be54d1adc917b55394e448f6fe`.
- TargetScanHuman release 8 coordinates are treated as transcript/UTR-relative
  and require an explicitly versioned annotation/assembly mapping before a
  genomic BED intersection.
- `consensus_hyb.py`, `extract_targeted_umi.py`,
  `targetscan_sites_to_bed12.py`, and the focused tests are
  optimization-authored MIT-covered files in this package. The targeted UMI
  script preserves the pinned route's R2-prefix transformation but requires an
  explicit integer length and records the declaring protocol because the
  upstream 9-nt prose and 10-nt default conflict.
