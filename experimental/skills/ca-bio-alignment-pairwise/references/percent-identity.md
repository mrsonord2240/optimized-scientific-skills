# Percent Identity Definitions (reference)

## Percent Identity: Definitions Matter

There are four common ways to calculate percent identity from the same alignment, producing different values:

| Method | Denominator | Best For |
|--------|-------------|----------|
| PID1 | Aligned positions + internal gaps | Gap-aware, conservative |
| PID2 | Aligned residue pairs (excluding gaps) | Always highest value |
| PID3 | Shorter sequence length | Length-normalized |
| PID4 | Mean sequence length | Best correlation with structural similarity |

**Practical impact**: Up to 11.5% difference between methods on a single alignment. Combined with different alignment algorithms, variation reaches 22%. Always report which method was used. The `counts()` method (SKILL.md, Alignment Counts) uses aligned non-gap positions (similar to PID2).
