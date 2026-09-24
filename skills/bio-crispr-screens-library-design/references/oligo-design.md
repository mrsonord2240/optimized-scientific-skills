# Oligo Design for Pooled Synthesis

## Oligo Design for Pooled Synthesis

**Goal:** Generate the final oligo sequence ready for chip-based synthesis. Vendor limits differ: Twist oligo pools cap at ~300 nt per oligo with no fixed pool size, GenScript's 92K format spans 20-170 nt, and Agilent OLS 244K spans 30-230 nt.

**Approach:** Add subpool PCR primers (so multiple sublibraries can share a synthesis array), the BsmBI/Esp3I overhang for golden-gate cloning into LentiGuide-Puro (Addgene 52963) or LentiCRISPRv2, and append the tracrRNA scaffold if the array length permits.

`scripts/build_oligo.py` builds the synthesis oligo (subpool forward primer + spacer + scaffold stub) and raises if it exceeds the 200 nt design budget:

```bash
python scripts/build_oligo.py GGATGGAGACGCATGATTCA --subpool 1
```

LentiGuide-Puro / LentiCRISPRv2 use BsmBI (Esp3I); the annealed-oligo overhangs are forward `5'-CACCG[spacer]-3'` and reverse `5'-AAAC[revcomp(spacer)]C-3'`. The scaffold stub is the first 33 nt of the Chen 2013 sgRNA(F+E) scaffold; lentiGuide-Puro (#52963) and lentiCRISPRv2 (#52961) carry the ORIGINAL scaffold, F+E belongs to lentiCRISPRv2-Opti (#163126).

**Subpool design:** A large synthesis pool can be partitioned into multiple sublibraries via subpool primers; each sub-PCR amplifies its subpool, allowing one synthesis batch to serve several screens. Typical subpool size: 10k-20k oligos.
