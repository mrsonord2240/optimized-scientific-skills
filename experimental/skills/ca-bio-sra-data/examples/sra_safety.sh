#!/bin/bash
# Shared safety helpers for public SRR/ERR/DRR download scripts.

require_sra_run_accession() {
    local accession="$1"
    if [[ ! "${accession}" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]; then
        echo "Invalid run accession: ${accession}. Expected uppercase SRR, ERR, or DRR followed only by digits." >&2
        return 2
    fi
}

make_owned_stage() {
    local out_dir="$1"
    local accession="$2"
    local purpose="$3"
    local out_root
    local stage

    require_sra_run_accession "${accession}" || return
    mkdir -p -- "${out_dir}"
    out_root="$(cd -- "${out_dir}" && pwd -P)" || return
    stage="$(mktemp -d "${out_root}/.${accession}.${purpose}.XXXXXX")" || return
    printf '%s\n' "${out_root}" > "${stage}/.sra-owned-stage"
    printf '%s\n' "${stage}"
}

cleanup_owned_stage() {
    local stage="$1"
    local out_dir="$2"
    local out_root

    out_root="$(cd -- "${out_dir}" && pwd -P)" || return 1
    case "${stage}" in
        "${out_root}"/.*) ;;
        *) echo "Refusing to clean a stage outside the requested output directory: ${stage}" >&2; return 1 ;;
    esac
    if [ ! -f "${stage}/.sra-owned-stage" ]; then
        echo "Refusing to clean an unmarked stage: ${stage}" >&2
        return 1
    fi
    rm -rf -- "${stage}"
}

publish_no_clobber() {
    local source_file="$1"
    local destination="$2"

    # -L catches a dangling symlink, for which -e alone is false. The hard-link
    # creation is the publication operation: it fails atomically if the name exists.
    if [ -e "${destination}" ] || [ -L "${destination}" ]; then
        echo "Refusing to overwrite existing output: ${destination}" >&2
        return 1
    fi
    if ! ln -- "${source_file}" "${destination}"; then
        echo "Refusing to publish over a concurrently created output: ${destination}" >&2
        return 1
    fi
    rm -f -- "${source_file}"
}
