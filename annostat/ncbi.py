"""Optional NCBI Datasets connectivity for accession-based inputs."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path


_ASSEMBLY_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$")


@dataclass(frozen=True, slots=True)
class FetchedGenome:
    """Describe one assembly extracted from an NCBI genome data package."""

    accession: str
    fasta: Path
    gff: Path
    metadata: dict[str, object]


def fetch_genomes(accessions: list[str], output_dir: Path) -> list[FetchedGenome]:
    """Download and extract FASTA/GFF3 pairs using NCBI's official CLI.

    The external ``datasets`` command owns networking, API-key handling, and ZIP
    validation. AnnStat only validates assembly accessions and exposes the two
    files needed by its local analysis workflows.
    """

    if not accessions:
        raise ValueError("at least one NCBI assembly accession is required")
    invalid = [accession for accession in accessions if not _ASSEMBLY_ACCESSION.fullmatch(accession)]
    if invalid:
        raise ValueError(
            "invalid NCBI assembly accession(s): " + ", ".join(invalid)
            + "; expected a versioned GCF_ or GCA_ accession"
        )
    executable = shutil.which("datasets")
    if executable is None:
        raise OSError(
            "NCBI Datasets CLI was not found; install it from "
            "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="annostat-ncbi-") as temporary_directory:
        archive = Path(temporary_directory) / "ncbi_dataset.zip"
        command = [
            executable, "download", "genome", "accession", *accessions,
            "--include", "genome,gff3", "--filename", str(archive), "--no-progressbar",
        ]
        result = subprocess.run(command, check=False, capture_output=True, text=True)
        if result.returncode:
            details = result.stderr.strip() or result.stdout.strip() or "unknown NCBI Datasets error"
            raise OSError(f"NCBI download failed: {details}")
        with zipfile.ZipFile(archive) as package:
            _safe_extract(package, output_dir)

    metadata_by_accession = _read_assembly_metadata(
        output_dir / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    )
    fetched = []
    for accession in accessions:
        directory = output_dir / "ncbi_dataset" / "data" / accession
        fasta_files = sorted(directory.glob("*_genomic.fna"))
        gff_files = sorted(directory.glob("*.gff")) + sorted(directory.glob("*.gff3"))
        if len(fasta_files) != 1 or len(gff_files) != 1:
            raise OSError(f"NCBI package for {accession} did not contain one genome FASTA and one GFF3")
        fetched.append(
            FetchedGenome(
                accession,
                fasta_files[0],
                gff_files[0],
                metadata_by_accession.get(accession, {"assembly_accession": accession}),
            )
        )
    return fetched


def _read_assembly_metadata(path: Path) -> dict[str, dict[str, object]]:
    """Select comparison-relevant fields from an NCBI assembly data report."""

    if not path.is_file():
        return {}
    metadata = {}
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError as error:
            raise OSError(f"invalid NCBI assembly metadata at line {line_number}: {error.msg}") from error
        accession = record.get("accession") or record.get("currentAccession")
        if not accession:
            continue
        organism = record.get("organism", {})
        assembly = record.get("assemblyInfo", {})
        annotation = record.get("annotationInfo", {})
        checkm = record.get("checkmInfo", {})
        metadata[str(accession)] = {
            "input_source": "NCBI Datasets",
            "assembly_accession": accession,
            "organism_name": organism.get("organismName"),
            "taxonomy_id": organism.get("taxId"),
            "assembly_name": assembly.get("assemblyName"),
            "assembly_level": assembly.get("assemblyLevel"),
            "annotation_provider": annotation.get("provider"),
            "annotation_release": annotation.get("releaseDate"),
            "annotation_pipeline": annotation.get("pipeline") or annotation.get("method"),
            "checkm_completeness": checkm.get("completeness"),
            "checkm_contamination": checkm.get("contamination"),
        }
    return metadata


def _safe_extract(package: zipfile.ZipFile, destination: Path) -> None:
    """Extract a ZIP only when every member remains below the destination."""

    destination = destination.resolve()
    for member in package.infolist():
        target = (destination / member.filename).resolve()
        if destination != target and destination not in target.parents:
            raise OSError(f"unsafe path in NCBI data package: {member.filename}")
    package.extractall(destination)
