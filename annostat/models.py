"""Data models shared by annostat modules."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Feature:
    """One feature line from a GFF3 file."""

    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float | None
    strand: str
    phase: int | None
    attributes: dict[str, str]

    @property
    def id(self) -> str:
        """Return a stable display identifier when GFF3 ID is absent."""

        return self.attributes.get(
            "ID", f"{self.seqid}:{self.start}-{self.end}:{self.strand}"
        )

    @property
    def length(self) -> int:
        """Return feature length for 1-based inclusive coordinates."""

        return self.end - self.start + 1


@dataclass(frozen=True, slots=True)
class CdsSequence:
    """Nucleotide and translated sequence for one CDS feature."""

    feature: Feature
    nucleotide: str
    coding_nucleotide: str
    protein: str
