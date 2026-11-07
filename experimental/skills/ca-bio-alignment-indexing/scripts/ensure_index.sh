#!/usr/bin/env bash
set -euo pipefail

ensure_index() {   # ensure_index file.bam|file.cram [extra samtools-index options, e.g. -@ 4]
    local f=$1; shift
    local stem=${f%.*} idx have=0 stale=0 csi=0 min_shift=14
    local -a candidates
    case $f in
        *.bam)  candidates=("$f.csi" "$stem.csi" "$f.bai" "$stem.bai") ;;
        *.cram) candidates=("$f.crai" "$stem.crai") ;;
        *) printf 'Expected a .bam or .cram file: %s\n' "$f" >&2; return 2 ;;
    esac
    for idx in "${candidates[@]}"; do
        [ -e "$idx" ] || continue
        have=1
        [ "$f" -nt "$idx" ] && stale=1
        case $idx in
            *.csi) csi=1
                    # CSI is BGZF-compressed; decompressed bytes 5-8 store min_shift.
                    min_shift=$(bgzip -cd "$idx" | od -An -j4 -N4 -tu4 | tr -d '[:space:]')
                    [ -n "$min_shift" ] || min_shift=14 ;;
        esac
    done
    if [ $have = 0 ] || [ $stale = 1 ]; then
        rm -f "${candidates[@]}"
        if [ $csi = 1 ]; then
            samtools index -c -m "$min_shift" "$@" "$f"  # preserve an existing CSI bin size
        else
            samtools index "$@" "$f"
        fi
    fi
}


for f in "$@"; do ensure_index "$f"; done
