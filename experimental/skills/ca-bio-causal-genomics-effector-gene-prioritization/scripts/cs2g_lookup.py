#!/usr/bin/env python3
"""
Reference: cS2G pre-computed SNP-to-gene scores (Gazal 2022 Nat Genet 54:827), Zenodo record
7754032, public and unauthenticated. Checked 2026-09-21 against cS2G_1000GEUR.zip (hg19,
9,997,231 SNPs, per-chromosome files `cS2G.<chr>.SGscore.gz`, columns SNP GENE cS2G INFO).
Stdlib only.

Look up cS2G gene scores for a list of rsIDs and rank the genes. cS2G links each SNP to genes
with scores that sum to <= 1 per SNP; INFO lists the constituent strategies that fired
(Promoter, ABC, EpiMap, Roadmap, GTeX_Finemapped, eQTLGen_Finemapped, ...).

Usage:
    python cs2g_lookup.py <cS2G_1000GEUR.zip | extracted dir> <chr> <rsid> [<rsid> ...]
    python cs2g_lookup.py cS2G_1000GEUR.zip 1 rs11206509 rs10788994

Get the data once (95 MB; download the zip, no install):
    curl -L -o cS2G_1000GEUR.zip https://zenodo.org/api/records/7754032/files/cS2G_1000GEUR.zip/content

Output: per-SNP rows (SNP, GENE, cS2G, INFO), then a per-locus table summing cS2G per gene
across the queried SNPs (a locus-level gene ranking, weighted by how many credible variants
link to each gene). Pass the credible-set variants, not the whole locus. cS2G is a per-SNP
aggregator; L2G is per-(locus, gene). Both are informative when they disagree.
"""
import gzip
import io
import os
import sys
import zipfile
from collections import defaultdict


def open_chr(source, chrom):
    name = f"cS2G.{chrom}.SGscore.gz"
    if os.path.isdir(source):
        for root, _, files in os.walk(source):
            if name in files:
                return gzip.open(os.path.join(root, name), "rt")
        raise FileNotFoundError(f"{name} not found under {source}")
    zf = zipfile.ZipFile(source)
    member = next((n for n in zf.namelist() if n.endswith("/" + name) and "__MACOSX" not in n), None)
    if member is None:
        raise FileNotFoundError(f"{name} not in {source}")
    return io.TextIOWrapper(gzip.open(zf.open(member), "rb"), encoding="utf-8")


def main(argv):
    if len(argv) < 4:
        sys.exit(__doc__)
    source, chrom, rsids = argv[1], argv[2], set(argv[3:])
    hits = []
    with open_chr(source, chrom) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        assert header[:3] == ["SNP", "GENE", "cS2G"], f"unexpected header: {header}"
        for line in fh:
            snp, gene, score, *info = line.rstrip("\n").split("\t")
            if snp in rsids:
                hits.append((snp, gene, float(score), info[0] if info else ""))
    missing = rsids - {h[0] for h in hits}
    print("SNP\tGENE\tcS2G\tINFO")
    for snp, gene, score, info in sorted(hits, key=lambda h: (h[0], -h[2])):
        print(f"{snp}\t{gene}\t{score:g}\t{info}")
    if missing:
        print(f"# not linked to any gene (or not in the 1000G EUR MAC>=5 panel): {sorted(missing)}",
              file=sys.stderr)
    per_gene = defaultdict(float)
    for _, gene, score, _ in hits:
        per_gene[gene] += score
    print("\n# locus-level ranking: sum of cS2G over queried SNPs")
    print("GENE\tsum_cS2G")
    for gene, total in sorted(per_gene.items(), key=lambda kv: -kv[1]):
        print(f"{gene}\t{total:g}")


if __name__ == "__main__":
    main(sys.argv)
