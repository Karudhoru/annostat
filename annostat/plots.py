"""Small dependency-free SVG chart writers."""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from html import escape
from pathlib import Path


def write_bar_chart(
    path: Path,
    title: str,
    counts: Mapping[str, int],
    *,
    sort_by_value: bool = True,
) -> None:
    """Write a labeled horizontal bar chart as SVG."""

    items = list(counts.items())
    if sort_by_value:
        items.sort(key=lambda item: (-item[1], item[0]))
    width = 900
    left = 130
    top = 70
    row_height = 34
    plot_width = width - left - 80
    height = max(180, top + row_height * len(items) + 45)
    maximum = max((value for _, value in items), default=1)
    elements = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width / 2}" y="35" text-anchor="middle" font-family="sans-serif" font-size="22" font-weight="bold">{escape(title)}</text>',
    ]
    for index, (label, value) in enumerate(items):
        y = top + index * row_height
        bar_width = plot_width * value / maximum
        elements.extend(
            [
                f'<text x="{left - 10}" y="{y + 20}" text-anchor="end" font-family="sans-serif" font-size="14">{escape(label)}</text>',
                f'<rect x="{left}" y="{y}" width="{bar_width:.1f}" height="24" rx="3" fill="#377eb8"/>',
                f'<text x="{left + bar_width + 8:.1f}" y="{y + 18}" font-family="sans-serif" font-size="13">{value}</text>',
            ]
        )
    elements.append("</svg>")
    path.write_text("\n".join(elements) + "\n", encoding="utf-8")


def write_histogram(path: Path, title: str, values: Sequence[int], bins: int = 20) -> None:
    """Write a basic SVG histogram for integer observations."""

    if not values:
        write_bar_chart(path, title, {})
        return
    minimum = min(values)
    maximum = max(values)
    span = max(1, maximum - minimum + 1)
    bin_width = max(1, (span + bins - 1) // bins)
    histogram: Counter[int] = Counter((value - minimum) // bin_width for value in values)
    labels = {
        f"{minimum + index * bin_width}-{min(maximum, minimum + (index + 1) * bin_width - 1)}": histogram[index]
        for index in range((span + bin_width - 1) // bin_width)
    }
    write_bar_chart(path, title, labels, sort_by_value=False)
