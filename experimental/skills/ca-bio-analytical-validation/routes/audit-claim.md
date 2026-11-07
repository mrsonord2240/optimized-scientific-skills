# Auditing a sensitivity claim

No script: work down the checklist and answer in words. For the numbers, call the other routes.

1. Is the VAF paired with input mass in GE? If not, compute lambda with `routes/ge-poisson.md`. Lambda below 3 means the template is often absent regardless of sequencing: 0.1% on 100 GE is noise.
2. Is the LoD tied to unique-molecule (consensus) depth and a replicate detection rate, from full-workflow replicates? A sampling calculation or k-of-N bound quoted as panel LoD95 is a theoretical bound (`routes/panel-integrated-lod.md`).
3. Does the error floor allow it? The floor is the per-base background rate, independent of sampling; the achieved LoD is the worse of the two ceilings.

| Error suppression | Typical per-base floor |
|-------------------|------------------------|
| Raw NGS | about 1e-3 |
| Single-strand UMI consensus | 1e-4 to 1e-5 |
| Duplex sequencing | below 1e-7 (theory below 1 per 1e9 nt) |

   Single-strand consensus does not remove template-resident damage (C to T deamination, G to T 8-oxoG): every copy of that strand inherits the lesion. Only duplex strand concordance catches it, so sub-1e-5 claims need duplex.
4. Is the background modelled per locus? Error depends on position and context; a flat global VAF cutoff loses sensitivity and specificity (integrated digital error suppression models it per locus and gives about 3-15x improvement, ctDNA to about 4e-5).
5. Detection or quantitation? Yes/no calls are governed by LoD; a reported VAF or tumour fraction needs LoQ (CV at or below 20%), and LoQ is never below LoD. Copy-number tumour-fraction estimation (ichorCNA, ULP-WGS) has a floor near 3%, which is a tumour-fraction limit, not a VAF.
6. Does the dilution-series LoD95 come with an interval and a bracketing design (`routes/lod95-probit.md`)?

Done when each of the six points is answered as supported, unsupported or not stated, and the claim is restated in the conditioned form "LoD95 at X ng, Y consensus depth, Z of N replicates".
