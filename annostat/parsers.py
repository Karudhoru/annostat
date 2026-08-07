"""Parsers for GFF3 annotations and FASTA sequences."""

from __future__ import annotations

from collections.abc import Iterator
from math import isfinite
from pathlib import Path
from urllib.parse import unquote

from annostat.models import Feature


def parse_attributes(text: str) -> dict[str, str]:
    """Parse the semicolon-separated ninth GFF3 column into a mapping.

    Percent-encoded keys and values are decoded according to GFF3 conventions.
    Repeated keys are retained as one comma-separated value rather than silently
    discarding earlier data.
    """

    attributes: dict[str, str] = {}
    if text == ".":
        return attributes

    for item in text.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if not separator:
            attributes[unquote(key)] = ""
            continue
        key = unquote(key)
        value = unquote(value)
        if key in attributes:
            # Preserve repeated attributes using GFF3's list separator.
            attributes[key] = f"{attributes[key]},{value}"
        else:
            attributes[key] = value
    return attributes


def parse_gff(path: str | Path) -> Iterator[Feature]:
    """Yield validated :class:`Feature` objects from a GFF3 file.

    Comment and directive lines are ignored, and parsing stops at an optional
    embedded FASTA section. Malformed field counts, coordinates, strands,
    scores, or phases raise ``ValueError`` with the source line number.
    """

    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if line == "##FASTA":
                # Sequence lines after this directive are not GFF3 features.
                break
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 tab-separated GFF3 fields, "
                    f"found {len(fields)}"
                )
            (
                seqid, source, feature_type, start_text, end_text, score_text,
                strand, phase_text, attribute_text,
            ) = fields
            if not seqid or not feature_type:
                raise ValueError(
                    f"{path}:{line_number}: seqid and feature type must not be empty"
                )
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: start and end must be integers"
                ) from error
            if start < 1 or end < start:
                raise ValueError(
                    f"{path}:{line_number}: invalid coordinate range {start}-{end}"
                )
            if strand not in {"+", "-", ".", "?"}:
                raise ValueError(f"{path}:{line_number}: invalid strand {strand!r}")

            try:
                score = None if score_text == "." else float(score_text)
                phase = None if phase_text == "." else int(phase_text)
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid score or phase"
                ) from error
            if score is not None and not isfinite(score):
                raise ValueError(f"{path}:{line_number}: score must be finite or .")
            if phase not in {None, 0, 1, 2}:
                raise ValueError(f"{path}:{line_number}: phase must be 0, 1, 2, or .")

            yield Feature(
                seqid=seqid,
                source=source,
                type=feature_type,
                start=start,
                end=end,
                score=score,
                strand=strand,
                phase=phase,
                attributes=parse_attributes(attribute_text),
            )


def parse_fasta(path: str | Path) -> dict[str, str]:
    """Read FASTA records into an uppercase sequence mapping.

    The first whitespace-delimited header token is used as the sequence ID so
    it matches the GFF3 ``seqid`` field. Empty and duplicate records are rejected
    with descriptive ``ValueError`` messages.
    """

    sequences: dict[str, list[str]] = {}
    current_id: str | None = None
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                current_id = header.split(maxsplit=1)[0]
                if current_id in sequences:
                    raise ValueError(
                        f"{path}:{line_number}: duplicate FASTA identifier {current_id!r}"
                    )
                sequences[current_id] = []
            elif current_id is None:
                raise ValueError(
                    f"{path}:{line_number}: sequence data appears before a FASTA header"
                )
            else:
                sequences[current_id].append("".join(line.split()).upper())

    if not sequences:
        raise ValueError(f"{path}: no FASTA records found")
    empty_identifiers = [identifier for identifier, parts in sequences.items() if not parts]
    if empty_identifiers:
        raise ValueError(
            f"{path}: FASTA record {empty_identifiers[0]!r} contains no sequence"
        )
    return {identifier: "".join(parts) for identifier, parts in sequences.items()}
