# Index formats

BAI covers coordinates below 2^29, CSI extends reach with depth sized from reference lengths, and CRAI indexes CRAM. CSI min_shift changes bin size rather than being the recommended way to increase reach. htslib prefers BAM CSI over BAI; stale competing indices may override fresh ones. For CRAM, REF_PATH is an MD5 reference cache rather than a plain FASTA directory.
