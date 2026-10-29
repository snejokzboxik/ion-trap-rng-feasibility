"""Fail-closed secret scan and reproducible manifest for the candidate public tree."""
from __future__ import annotations
import hashlib, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE_NAMES = {"LITERATURE_SESSION_NOTES.json", "PUBLICATION_MANIFEST.json"}
EXCLUDE_PARTS = {".git", ".venv", "venv", "__pycache__", ".pytest_cache", "build", "dist"}
SECRET_PATTERNS = {
    "private_key": re.compile(rb"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----"),
    "github_token": re.compile(rb"\bgh[pousr]_[A-Za-z0-9_]{20,}\b"),
    "aws_access_key": re.compile(rb"\bAKIA[0-9A-Z]{16}\b"),
    "generic_secret_assignment": re.compile(rb"(?i)(api[_-]?key|access[_-]?token|client[_-]?secret|password)\s*[:=]\s*['\"][^'\"\s]{12,}"),
    "private_windows_path": re.compile(rb"(?i)[A-Z]:\\Users\\[A-Za-z0-9._-]+\\"),
}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

files, findings = [], []
for path in sorted(ROOT.rglob("*")):
    if not path.is_file():
        continue
    rel = path.relative_to(ROOT)
    if (path.name in EXCLUDE_NAMES or any(part in EXCLUDE_PARTS for part in rel.parts)
            or any(part.startswith("pytest_") for part in rel.parts)):
        continue
    data = path.read_bytes()
    for label, pattern in SECRET_PATTERNS.items():
        if pattern.search(data):
            findings.append({"path": rel.as_posix(), "rule": label})
    files.append({"path": rel.as_posix(), "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
if findings:
    raise SystemExit("Publication audit failed: " + json.dumps(findings, ensure_ascii=False))
excluded = []
for name, reason in [("LITERATURE_SESSION_NOTES.json", "raw search-response snippets and third-party excerpts; redistribution rights not assessed; curated audit is retained"),
                     ("PUBLICATION_MANIFEST.json", "manifest does not hash itself")]:
    p = ROOT / name
    if p.exists():
        excluded.append({"path": name, "size_bytes": p.stat().st_size, "sha256": sha256(p), "reason": reason})
manifest = {
    "project": "ion-trap-rng-feasibility",
    "status": "published public repository; verified fresh clone and tests during Stage R0.1",
    "publication": {
        "url": "https://github.com/snejokzboxik/ion-trap-rng-feasibility",
        "visibility": "PUBLIC",
        "initial_published_revision": "aa73cfb439f3f439bcbf2087ff983c53162c8b28",
    },
    "license": "No license file or license grant included; rights compatibility not established.",
    "files": files, "excluded": excluded,
}
(ROOT / "PUBLICATION_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
print(json.dumps({"included_files": len(files), "included_bytes": sum(x["size_bytes"] for x in files), "excluded": excluded, "secret_findings": 0}, indent=2, ensure_ascii=False))
