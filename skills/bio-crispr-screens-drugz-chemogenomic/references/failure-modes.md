# drugZ Failure Modes

### drugZ shows no synthetic-lethal hits despite known sensitizing genes

**Trigger:** Comparing drug vs Day-0 instead of drug vs vehicle.
**Mechanism:** Day-0 comparison conflates drug effect with normal-culture proliferation; essential genes drop in both conditions, masking drug-specific sensitization.
**Symptom:** PARPi screen shows no sensitization at BRCA1/BRCA2 despite expected biology.
**Fix:** Re-run with vehicle samples passed to `-c`. The drug-vs-vehicle is the canonical comparison.

### High false-positive rate among essential genes

**Trigger:** Essential genes drop out in both vehicle and drug arms; small relative shift gives misleadingly high Z.
**Mechanism:** drugZ's Z-score is symmetric; essential genes drop in both arms but slightly more in drug -> "synthetic lethal" call.
**Symptom:** Hit list dominated by RPS, RPL, EIF essentials.
**Fix:** Use `-r` with a comma-delimited list of essential gene names to exclude (see "Removing Reference Genes"), then check that the excluded genes really are missing from the output; or filter the output post-hoc.

### Unstable hits across libraries or sub-samples

**Trigger:** Insufficient sgRNAs per gene; small effect sizes.
**Mechanism:** drugZ's per-gene sumZ depends on enough sgRNAs to be stable; with 3-4 sgRNAs/gene, single-guide noise drives the ranking.
**Symptom:** Top hits move when you re-sequence, sub-sample replicates, or switch library.
**Note:** Re-running drugZ on the same input cannot show this. drugZ has no sampling step and no seed: identical input gives a byte-identical output file, so a rerun is not a stability check.
**Fix:** Use a 6+ sgRNAs/gene library (Avana, Dolcetto); check stability by holding out a replicate or bootstrapping the guides yourself; or use MAGeCK MLE.

### drugZ ignores dose information

**Trigger:** Multi-dose screen analyzed at highest dose only.
**Mechanism:** drugZ doesn't model dose; running at one dose loses the dose-response information.
**Symptom:** Hits at high dose may be dose-specific (not true responders).
**Fix:** Run drugZ at each dose; require consistency across doses for high-confidence hits.

### Drug-target gene appears as "suppressor"

**Trigger:** Loss of drug target reduces drug binding, increasing drug resistance.
**Mechanism:** Real biology -- drug target itself is a resistance gene from a KO perspective.
**Symptom:** Drug-target gene like PARP1 appears in suppressor list for PARPi screen.
**Fix:** Expected biology. Annotate the drug target separately. The suppressor list is correct.

### Low replicate concordance (Pearson < 0.85) or a crashed run

**Trigger:** replicate Pearson (within an arm, on log2 counts) below 0.85, or drugZ exits with a traceback.
**Mechanism:** drugZ's small-effect sensitivity also amplifies noise, so a noisy replicate produces false positives.
**Fix:** run `screen-qc` first; drop the worst replicate from `-c`/`-x` (drugZ accepts any number of columns, checked on 2 vs 2) and re-run; if hits change materially, flag the screen as QC-failed instead of reporting them. For a crash, check the guide count against `--half_window_size` and that every `-c`/`-x` name matches a header column exactly.
