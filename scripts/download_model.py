"""Download the exact Yeast9 model this project was built against.

The SBML file is ~12 MB, so it is not committed to the repository. Instead we
pin the release tag and the SHA-256 of the file, and verify the download. That
way a clean clone reproduces the analysis against *the same* model, and any
silent change upstream is caught immediately.

Usage:
    python scripts/download_model.py
"""

import hashlib
import json
import pathlib
import sys
import urllib.request

# Provenance is kept in a committed JSON file so it is version-controlled
# alongside the code and can be cited in the README.
PROVENANCE = pathlib.Path(__file__).resolve().parents[1] / "data" / "model" / "MODEL_PROVENANCE.json"


def sha256_of(path: pathlib.Path) -> str:
    """Return the SHA-256 hex digest of a file, read in chunks."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    meta = json.loads(PROVENANCE.read_text())
    target = PROVENANCE.parent / pathlib.Path(meta["file"]).name

    # Already present and correct? Then there is nothing to do.
    if target.exists() and sha256_of(target) == meta["sha256"]:
        print(f"{target} already present, checksum matches.")
        return 0

    print(f"Downloading Yeast9 {meta['release_tag']} from {meta['source_url']}")
    target.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(meta["source_url"], target)

    actual = sha256_of(target)
    if actual != meta["sha256"]:
        # Refuse to leave a file we cannot vouch for: the analysis would be
        # reproducible in name only.
        target.unlink(missing_ok=True)
        print(
            "CHECKSUM MISMATCH -- download discarded.\n"
            f"  expected {meta['sha256']}\n  got      {actual}",
            file=sys.stderr,
        )
        return 1

    size_mb = target.stat().st_size / 1e6
    print(f"OK: {target} ({size_mb:.1f} MB), checksum matches {meta['sha256'][:12]}...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
