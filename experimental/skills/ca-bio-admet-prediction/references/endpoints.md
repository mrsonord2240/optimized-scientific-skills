## hERG Cardiotoxicity Endpoint

hERG (KCNH2) blockade can contribute to QT prolongation and Torsades de Pointes and is an important non-clinical cardiac-safety endpoint. Follow the current ICH S7B/E14 and regulator-specific guidance applicable to the program rather than treating one model output as a regulatory conclusion.

| Model | Architecture | Training data | AUC | Reference |
|-------|--------------|---------------|-----|-----------|
| Shan et al. D-MPNN + MOE | D-MPNN + 206 MOE descriptors | 7,889 compounds | 0.956 (random split) | Shan 2022 |
| CardioTox-net | Five DL base representations + neural meta-ensemble | BindingDB, ChEMBL, and literature | 0.930 (10-fold meta-validation) | Karim 2021 |
| ADMETlab 3.0 hERG | DMPNN multi-task | Internal | 0.92 (reported) | Fu 2024 |
| ProTox-3.0 cardiotoxicity | RF-based classifier | 5,252 ChEMBL compounds with hERG IC50/Ki | 0.86 CV; 0.95 external | Banerjee 2024 |

**Interpretation:** A single-model probability > 0.5 is NOT a kill signal. Triangulate multiple hERG-specific models and a literature search. ProTox-3.0 calls the endpoint cardiotoxicity, but its model specifically predicts small-molecule hERG blockers; it should not be treated as an independent non-hERG mechanism. Consider exposure relative to measured hERG potency and confirm important decisions experimentally rather than applying a universal safe/unsafe IC50 cutoff.

## CYP Inhibition (DDI Risk)

5 CYP isoforms cover most clinically relevant DDIs:

| CYP | Substrates (drugs) | Inhibitor flag if predicted prob | Action |
|-----|--------------------|-----------------------------------|--------|
| CYP3A4 | many drug classes | Model-specific threshold | Interpret inhibitor and substrate assays separately |
| CYP2D6 | beta-blockers, antidepressants | Model-specific threshold | Include polymorphism and exposure context |
| CYP2C9 | warfarin, NSAIDs | Model-specific threshold | Evaluate clinical substrate/exposure context |
| CYP2C19 | PPIs, clopidogrel | Model-specific threshold | Include polymorphism and assay context |
| CYP1A2 | caffeine, theophylline | Model-specific threshold | Include induction, diet, and smoking context |
