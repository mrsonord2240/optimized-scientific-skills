"""Materialise the two Skill arms as plain directory trees, pinned and verified.

  arms/upstream/skills/<id>/   from GPTomics/bioSkills @ d91ed3d   (tarball in arms/_up)
  arms/optimized/skills/<id>/  from optimized-scientific-skills @ HEAD (must be clean)

Both use the published id as the directory name so the two arms are indistinguishable
by path. The 'none' arm has no skills dir at all.
"""
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
OPT = Path(r"F:\optimized-scientific-skills")
UP = next((HERE / "arms" / "_up").iterdir())
OUT = HERE / "arms"

prov = json.loads((OPT / "PROVENANCE.json").read_text(encoding="utf-8"))
opt_commit = subprocess.check_output(["git", "-C", str(OPT), "rev-parse", "HEAD"], text=True).strip()
dirty = subprocess.check_output(["git", "-C", str(OPT), "status", "--porcelain", "--", "skills", "PROVENANCE.json"], text=True).strip()
assert not dirty, f"optimized repo has uncommitted changes:\n{dirty}"
up_commit = prov["sources"]["gptomics-bioskills"]["upstream_commit"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def tree(d: Path):
    return {str(p.relative_to(d)).replace("\\", "/"): sha(p) for p in sorted(d.rglob("*")) if p.is_file()}


for arm in ("upstream", "optimized"):
    shutil.rmtree(OUT / arm, ignore_errors=True)

manifest = {"upstream_commit": up_commit, "optimized_commit": opt_commit, "skills": {}}
for s in prov["skills"]:
    sid = s["id"]
    src_up = UP / s["upstream_path"]
    src_opt = OPT / "skills" / sid
    assert src_up.is_dir(), f"missing upstream dir for {sid}: {src_up}"
    assert src_opt.is_dir(), f"missing optimized dir for {sid}"
    shutil.copytree(src_up, OUT / "upstream" / "skills" / sid)
    shutil.copytree(src_opt, OUT / "optimized" / "skills" / sid)
    a, b = tree(OUT / "upstream" / "skills" / sid), tree(OUT / "optimized" / "skills" / sid)
    manifest["skills"][sid] = {
        "relative_to_upstream": s["relative_to_upstream"],
        "domain": s["upstream_path"].split("/")[0],
        "files_identical": sorted(k for k in a if k in b and a[k] == b[k]),
        "files_changed": sorted(k for k in a if k in b and a[k] != b[k]),
        "files_only_upstream": sorted(k for k in a if k not in b),
        "files_only_optimized": sorted(k for k in b if k not in a),
    }

# Sanity: 'licence-declaration-only' skills must differ only by SKILL.md (the license line).
bad = []
for sid, m in manifest["skills"].items():
    if m["relative_to_upstream"] == "licence-declaration-only":
        if m["files_changed"] != ["SKILL.md"] or m["files_only_upstream"] or m["files_only_optimized"]:
            bad.append(sid)
        else:
            u = (OUT / "upstream" / "skills" / sid / "SKILL.md").read_text(encoding="utf-8").splitlines()
            o = (OUT / "optimized" / "skills" / sid / "SKILL.md").read_text(encoding="utf-8").splitlines()
            extra = [l for l in o if l not in u]
            gone = [l for l in u if l not in o]
            if extra != ["license: MIT"] and not (all(l.startswith("license:") for l in extra) and not gone):
                bad.append(f"{sid} (SKILL.md diff beyond license: {extra[:2]} / {gone[:2]})")
manifest["licence_only_verification_failures"] = bad
(OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")

n_mod = sum(1 for m in manifest["skills"].values() if m["relative_to_upstream"] == "modified")
print(f"upstream {up_commit[:7]}  optimized {opt_commit[:7]}  skills {len(manifest['skills'])}  modified {n_mod}")
print("licence-only verification failures:", bad or "none")
