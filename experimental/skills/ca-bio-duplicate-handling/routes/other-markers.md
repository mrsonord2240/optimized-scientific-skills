# Other duplicate markers

Commands for samblaster (streaming after `bwa-mem2`), Picard `MarkDuplicates`, biobambam2, sambamba, mapDamage `--rescale` for ancient DNA, and `pbmarkdup` are in `references/alternative-markers.md`. Read it, copy the command for the tool asked for, and run it. These have no wrapper script.

- Picard reports `READ_PAIR_DUPLICATES` in pairs; samtools counts reads (50 pairs = 100 reads).
- Picard needs an `@RG` line (`bwa-mem2 mem -R`, or `samtools addreplacerg`) or it throws a NullPointerException.
- Ancient DNA: mark first, then mapDamage `--rescale` before variant calling.

Done when the marked BAM exists and the tool's duplicate count is reported with its units.
