"""Optional filters for exporting selected coding sequences."""

from __future__ import annotations

from dataclasses import dataclass

from annostat.analysis import cog_categories
from annostat.models import Feature


@dataclass(frozen=True, slots=True)
class CdsFilter:
    """Describe optional criteria applied to filtered CDS exports.

    Filtering never changes the complete analysis or its standard output files.
    When at least one criterion is active, matching records are additionally
    written beneath the ``filtered`` output directory.
    """

    min_length: int | None = None
    max_length: int | None = None
    require_cog: bool = False
    exclude_hypothetical: bool = False

    @property
    def active(self) -> bool:
        """Return whether at least one filtering criterion is enabled."""

        return any(
            (
                self.min_length is not None,
                self.max_length is not None,
                self.require_cog,
                self.exclude_hypothetical,
            )
        )

    def matches(self, feature: Feature, *, length: int | None = None) -> bool:
        """Return whether a CDS satisfies every enabled criterion.

        ``length`` lets multipart callers supply the joined CDS length while
        preserving the original single-feature API.
        """

        if feature.type != "CDS":
            return False
        observed_length = feature.length if length is None else length
        if self.min_length is not None and observed_length < self.min_length:
            return False
        if self.max_length is not None and observed_length > self.max_length:
            return False
        if self.require_cog and not cog_categories(feature):
            return False
        if self.exclude_hypothetical and "hypothetical protein" in feature.attributes.get(
            "product", ""
        ).lower():
            return False
        return True

    def as_dict(self) -> dict[str, int | bool | None]:
        """Return a JSON-safe representation for reports and provenance."""

        return {
            "active": self.active,
            "min_length": self.min_length,
            "max_length": self.max_length,
            "require_cog": self.require_cog,
            "exclude_hypothetical": self.exclude_hypothetical,
        }
