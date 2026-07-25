"""Download and verify the public cBioPortal TCGA-LAML study archive."""

from __future__ import annotations

import hashlib
import shutil
import tarfile
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
ARCHIVE = RAW / "laml_tcga_pub.tar.gz"
EXTRACTED = RAW / "laml_tcga_pub"
URL = "https://datahub.assets.cbioportal.org/laml_tcga_pub.tar.gz"
EXPECTED_SHA256 = "be0b4bb8a0481acd6e4127f4e0a15889de9eb7fd4bc6e663cc90f65db419b0b9"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_extract(archive: tarfile.TarFile, destination: Path) -> None:
    destination_resolved = destination.resolve()
    for member in archive.getmembers():
        target = (destination / member.name).resolve()
        if destination_resolved not in target.parents and target != destination_resolved:
            raise RuntimeError(f"Unsafe archive member: {member.name}")
    archive.extractall(destination, filter="data")


def main() -> None:
    RAW.mkdir(parents=True, exist_ok=True)

    if not ARCHIVE.exists():
        partial = ARCHIVE.with_suffix(ARCHIVE.suffix + ".part")
        print(f"Downloading {URL}")
        with urllib.request.urlopen(URL) as response, partial.open("wb") as output:
            shutil.copyfileobj(response, output)
        partial.replace(ARCHIVE)

    observed = sha256(ARCHIVE)
    if observed.lower() != EXPECTED_SHA256:
        raise RuntimeError(
            "Source archive hash differs from the recorded analysis version.\n"
            f"Expected: {EXPECTED_SHA256}\n"
            f"Observed: {observed}\n"
            "Do not continue without reviewing the upstream version change."
        )

    print(f"Verified SHA-256: {observed}")
    if EXTRACTED.exists():
        print(f"Already extracted: {EXTRACTED}")
        return

    with tarfile.open(ARCHIVE, "r:gz") as archive:
        safe_extract(archive, RAW)
    print(f"Extracted to: {EXTRACTED}")


if __name__ == "__main__":
    main()

