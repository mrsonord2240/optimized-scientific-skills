#!/bin/bash
# Shared fallback for public runs when the SRA Toolkit path cannot complete.
# Source this file from a script in this directory; do not run it directly.

run_ena_fallback() {
    local accession="$1"
    local out_dir="$2"
    local script_dir="$3"
    local stage_dir
    local accessions_file
    local files
    local destination

    # shellcheck source=sra_safety.sh
    source "${script_dir}/sra_safety.sh"
    require_sra_run_accession "${accession}" || return
    stage_dir=$(make_owned_stage "${out_dir}" "${accession}" "ena-fallback") || return
    accessions_file="${stage_dir}/accessions.txt"
    printf '%s\n' "${accession}" > "${accessions_file}"

    echo "=== SRA Toolkit route could not complete; trying verified ENA mirror fallback ==="
    if ! bash "${script_dir}/download_batch.sh" "${accessions_file}" "${stage_dir}"; then
        cleanup_owned_stage "${stage_dir}" "${out_dir}"
        return 1
    fi
    if [ -s "${stage_dir}/failed.txt" ] || ! compgen -G "${stage_dir}/*.fastq.gz" >/dev/null; then
        echo "ENA fallback did not produce verified FASTQ; no files were published." >&2
        if [ -s "${stage_dir}/failed.txt" ]; then
            echo "Failed accession(s): $(tr '\n' ' ' < "${stage_dir}/failed.txt")" >&2
        fi
        cleanup_owned_stage "${stage_dir}" "${out_dir}"
        return 1
    fi

    shopt -s nullglob
    files=("${stage_dir}"/*.fastq.gz)
    shopt -u nullglob
    for source_file in "${files[@]}"; do
        destination="${out_dir}/$(basename "${source_file}")"
        if [ -e "${destination}" ] || [ -L "${destination}" ]; then
            echo "Refusing to overwrite existing output: ${destination}" >&2
            cleanup_owned_stage "${stage_dir}" "${out_dir}"
            return 1
        fi
    done
    for source_file in "${files[@]}"; do
        destination="${out_dir}/$(basename "${source_file}")"
        if ! publish_no_clobber "${source_file}" "${destination}"; then
            cleanup_owned_stage "${stage_dir}" "${out_dir}"
            return 1
        fi
    done
    cleanup_owned_stage "${stage_dir}" "${out_dir}"
}
