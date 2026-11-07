# Other logo encodings

```bash
weblogo --format pdf --sequence-type dna --color-scheme classic --units bits --composition equiprobable --fineprint '' --size large < aligned.fa > logo.pdf
```

Use the matching sequence type, including RNA when applicable, and document WebLogo's prior or `--weight` choice. The source `examples/seqlogo_phd.R` demonstrates relative-entropy and uniform-background logos, multi-logo comparisons, protein coloring and PWM input using bundled fixtures; run from `examples/` and replace inputs before interpreting. It writes fixed output names. For Logomaker, positions are rows and letters columns; explicit background and pseudocounts matter, and signed weight logos are not an EDLogo implementation. Protein, signed and arbitrary alphabet workflows need adaptation of the supplied examples or tool APIs; there is no generic CLI wrapper for them. Compare different alphabets on a common normalized scale or separately.

Done when input orientation, alphabet, method, prior and actual output are verified.
