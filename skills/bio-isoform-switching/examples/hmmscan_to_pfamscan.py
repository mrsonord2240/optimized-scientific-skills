#!/usr/bin/env python3
"""Convert `hmmscan --domtblout` output to the pfam_scan.pl-style table that
IsoformSwitchAnalyzeR::analyzePFAM() reads.

analyzePFAM() rejects a raw --domtblout file ("more columns than column names") and recognises a
pfam_scan file only if some row carries a clan accession (CL...) in column 15, so the clan of each
Pfam family is added from Pfam-A.clans.tsv(.gz).

    hmmscan --cut_ga --domtblout pfam_domtbl.txt Pfam-A.hmm isoformSwitchAnalyzeR_isoform_AA.fasta > /dev/null
    curl -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
    python hmmscan_to_pfamscan.py pfam_domtbl.txt pfam_scanfmt.txt Pfam-A.clans.tsv.gz

Checked with HMMER 3.4 domtblout (Pfam-A release of 2026-09) and IsoformSwitchAnalyzeR 2.6.0.
Run hmmscan with --cut_ga: every row written here is marked significant (column 14 = 1), which is only
true after the Pfam gathering threshold.
"""
import gzip
import io
import sys


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else io.open(path, encoding="utf-8")


def main(domtbl, out, clans_path):
    clan = {}
    with open_text(clans_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 1 and f[1]:
                clan[f[0]] = f[1]  # Pfam accession (no version) -> clan accession
    if not clan:
        sys.exit("no clan rows read from %s (expected Pfam-A.clans.tsv: accession<TAB>clan<TAB>...)" % clans_path)
    rows = []
    n_clan = 0
    with io.open(domtbl, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 21:
                sys.exit("not an hmmscan --domtblout line (%d fields): %r" % (len(f), line[:80]))
            hmm_name, hmm_acc, hmm_len, query = f[0], f[1], int(f[2]), f[3]
            i_evalue, score = f[12], f[13]
            hmm_from, hmm_to, ali_from, ali_to, env_from, env_to = f[15:21]
            cl = clan.get(hmm_acc.split(".")[0], "No_clan")
            n_clan += cl != "No_clan"
            rows.append("%s %6s %6s %6s %6s %s  %s  Domain %5s %5s %5d %8s %10s 1 %s  "
                        % (query, ali_from, ali_to, env_from, env_to, hmm_acc, hmm_name, hmm_from, hmm_to, hmm_len, score, i_evalue, cl))
    if not n_clan:
        sys.stderr.write("warning: no row has a clan; analyzePFAM() will refuse this file\n")
    with io.open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")
    print("wrote %d domain rows to %s" % (len(rows), out))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
