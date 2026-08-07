"""Data models shared by annostat modules."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Feature:
    """Represent one parsed feature line from a GFF3 annotation.

    Coordinates remain 1-based and inclusive, matching the GFF3 file. Optional
    score and phase values are represented by ``None`` when the source uses ``.``.
    """

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
    """Store the sequence forms derived from one CDS feature.

    ``nucleotide`` contains the complete strand-oriented feature sequence,
    whereas ``coding_nucleotide`` has the GFF3 phase removed and is the sequence
    used to produce ``protein``.
    """

    feature: Feature
    nucleotide: str
    coding_nucleotide: str
    protein: str
