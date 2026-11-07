# betaPrior deprecation timeline

## betaPrior Deprecation Timeline

DESeq2 v1.0-v1.15: `betaPrior = TRUE` was the default. Shrinkage was baked in and p-values were computed on the shrunken estimates. v1.16 (2017) flipped the default to `FALSE` and introduced `lfcShrink()`. Today: `betaPrior = TRUE` is quasi-deprecated and is the ONLY path to a p-value of the shrunken estimate. Most users do not want this. `lfcShrink()` with apeglm and the unshrunken Wald p-value is the modern compromise.

The contrast= vs name= numerical difference that older tutorials describe was specific to `betaPrior = TRUE`. With the current default, `name=` and `contrast=` for the same comparison return identical LFCs.
