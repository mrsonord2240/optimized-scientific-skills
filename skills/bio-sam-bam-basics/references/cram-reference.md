## CRAM Reference Resolution (Critical)

CRAM stores reads relative to a reference; without it, the file is unreadable. htslib resolves the reference in this order:

1. Command-line `-T ref.fa` / `--reference`
2. `REF_CACHE` env var (local MD5-named cache; searched *before* `REF_PATH`)
3. `REF_PATH` env var (colon-separated; each element matched by the `@SQ M5:` MD5). A remote server such as EBI ENA is consulted only if its URL is present here -- it was the built-in default through htslib 1.21, but that default was **removed in 1.22** to reduce EBI load, so modern htslib does no network lookup unless that URL is added explicitly.
4. Local file named in the `@SQ UR:` header tag (local / `file://` paths only; `http`/`ftp` URIs in `UR:` are ignored)

On HPC nodes without internet, populate a local cache once:
```bash
mkdir -p $HOME/cram_cache
seq_cache_populate.pl -root $HOME/cram_cache reference.fa
export REF_CACHE=$HOME/cram_cache/%2s/%2s/%s
export REF_PATH=$REF_CACHE   # local only; no network/ENA lookup

samtools quickcheck -v file.cram                 # header + EOF only: passes with the reference missing and with a corrupt body
samtools view -o /dev/null file.cram && echo ok  # full decode: exit 1 if the reference cannot be resolved or a slice is corrupt
```
`samtools view -c`, `flagstat` and `idxstats` never decode bases, so they succeed on a CRAM whose reference is unreachable (checked on 1.24); only a full decode (`view -o /dev/null`, or `stats`) proves it.

`samtools view -C` without `-T` does not fail when no reference can be resolved: it warns (`Enabling embed_ref=2`) and embeds the read sequences. For a self-contained CRAM that decodes with no reference, use `-T ref.fa --output-fmt-option embed_ref=1`.

CRAM can be made irreversibly lossy, but the `archive` profile is NOT how: `--output-fmt-option archive` is a *lossless* maximum-compression preset (fqzcomp quality codec, name tokenization, larger slices) that does not alter bases or qualities. Irreversible loss comes instead from explicit quality **binning** (e.g. Illumina 8-bin), which must be applied deliberately and is harmful for low-coverage / somatic / forensic / archival data. Convert against the *exact* reference the BAM was aligned to (matched by `@SQ M5:`). A wrong reference is refused, not silently accepted: each slice stores its reference MD5, so decoding fails with `MD5 checksum reference mismatch` (exit 1). Bases are silently wrong only with `--input-fmt-option ignore_md5=1`.
