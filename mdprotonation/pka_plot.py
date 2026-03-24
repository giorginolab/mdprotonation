from __future__ import annotations

from dataclasses import dataclass

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .app_state import site_state_sort_key
from .charge_colors import charge_color_band, charge_legend_bands
from .presentation import TRANSITION_MARKER, TRANSITION_WINDOW_LABEL
from .protonation import SiteState, site_charge_at_ph

PH_MIN = 0.0
PH_MAX = 14.0
LANE_COLOR = charge_color_band(0.0).background
LANE_EDGE_COLOR = "#d1d5db"
CURRENT_PH_COLOR = "#3f3f46"
TRANSITION_WINDOW_COLOR = "#6b7280"
PH_SEGMENT_STEP = 0.1


@dataclass(frozen=True)
class PkaPlotRow:
    label: str
    pka: float
    pka_marker: float
    charge_segments: tuple["PkaPlotSegment", ...]
    is_transitioning: bool


@dataclass(frozen=True)
class PkaPlotSegment:
    start_ph: float
    end_ph: float
    color: str


def build_pka_plot_rows(site_states: list[SiteState]) -> list[PkaPlotRow]:
    rows: list[PkaPlotRow] = []
    sorted_states = sorted(site_states, key=site_state_sort_key)
    for state in sorted_states:
        rows.append(
            PkaPlotRow(
                label=_format_plot_label(state),
                pka=state.site.pka,
                pka_marker=_clip_to_plot(state.site.pka, inset=0.22),
                charge_segments=_build_charge_segments(state),
                is_transitioning=state.dominant_state == "Transitioning",
            )
        )
    return rows


def create_pka_plot_figure(site_states: list[SiteState], current_ph: float) -> Figure:
    rows = build_pka_plot_rows(site_states)
    figure_height = max(4.5, len(rows) * 0.34 + 1.4)
    figure, axis = plt.subplots(figsize=(12.5, figure_height))

    for index, row in enumerate(rows):
        axis.barh(
            index,
            PH_MAX - PH_MIN,
            left=PH_MIN,
            height=0.82,
            color="none",
            edgecolor=LANE_EDGE_COLOR,
            linewidth=0.7,
            zorder=1,
        )
        for segment in row.charge_segments:
            axis.barh(
                index,
                segment.end_ph - segment.start_ph,
                left=segment.start_ph,
                height=0.82,
                color=segment.color,
                edgecolor=segment.color,
                linewidth=0.0,
                zorder=2,
            )
        axis.text(
            row.pka_marker,
            index,
            f"{row.pka:.2f}",
            ha="center",
            va="center",
            fontsize=8,
            color="#1f2937",
            bbox={
                "boxstyle": "round,pad=0.12",
                "facecolor": "#f8fafc",
                "edgecolor": "none",
                "alpha": 0.72,
            },
            zorder=4,
        )

    for boundary in (current_ph - 1.0, current_ph + 1.0):
        if PH_MIN <= boundary <= PH_MAX:
            axis.axvline(
                boundary,
                color=TRANSITION_WINDOW_COLOR,
                linewidth=1.6,
                alpha=0.55,
                linestyle="--",
                zorder=4,
            )
    axis.axvline(
        current_ph,
        color=CURRENT_PH_COLOR,
        linewidth=2.4,
        alpha=0.92,
        zorder=5,
    )
    axis.set_xlim(PH_MIN, PH_MAX)
    axis.set_ylim(-0.8, len(rows) - 0.2)
    axis.set_yticks(range(len(rows)))
    axis.set_yticklabels([row.label for row in rows], fontsize=9)
    axis.invert_yaxis()

    axis.xaxis.tick_top()
    axis.xaxis.set_label_position("top")
    axis.set_xlabel("pKa", fontsize=12)
    axis.set_xticks(range(0, 15, 2))
    axis.set_xticks(range(0, 15), minor=True)
    axis.grid(axis="x", which="major", color="#374151", alpha=0.24, linewidth=1.0)
    axis.grid(
        axis="x",
        which="minor",
        color="#6b7280",
        alpha=0.18,
        linestyle="--",
        linewidth=0.8,
    )
    axis.tick_params(axis="y", length=0)
    axis.tick_params(axis="x", labelsize=10)

    for spine_name in ("left", "right", "bottom"):
        axis.spines[spine_name].set_visible(False)
    axis.spines["top"].set_color("#6b7280")
    axis.spines["top"].set_linewidth(1.0)

    legend_handles = [
        Patch(facecolor=band.color, edgecolor=LANE_EDGE_COLOR, label=band.label)
        for band in charge_legend_bands()
    ] + [
        Line2D(
            [0],
            [0],
            color=CURRENT_PH_COLOR,
            linewidth=2.4,
            label=f"Selected pH {current_ph:.1f}",
        ),
        Line2D(
            [0],
            [0],
            color=TRANSITION_WINDOW_COLOR,
            linewidth=1.6,
            linestyle="--",
            label=TRANSITION_WINDOW_LABEL,
        ),
    ]
    axis.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=4,
        frameon=False,
        fontsize=9,
    )

    figure.tight_layout()
    return figure


def _format_plot_label(state: SiteState) -> str:
    site = state.site
    insertion_code = site.insertion_code.strip()
    residue_token = f"{site.residue_number}{insertion_code}" if insertion_code else str(
        site.residue_number
    )
    prefix = f"{TRANSITION_MARKER} " if state.dominant_state == "Transitioning" else ""
    return f"{prefix}{site.chain_id}:{residue_token}-{site.residue_type}"


def _build_charge_segments(state: SiteState) -> tuple[PkaPlotSegment, ...]:
    segments: list[PkaPlotSegment] = []
    current_start = PH_MIN
    current_color = charge_color_band(
        site_charge_at_ph(state.site, PH_MIN + (PH_SEGMENT_STEP / 2.0))
    ).background

    ph = PH_MIN
    while ph < PH_MAX:
        start_ph = ph
        end_ph = min(ph + PH_SEGMENT_STEP, PH_MAX)
        midpoint = start_ph + ((end_ph - start_ph) / 2.0)
        color = charge_color_band(site_charge_at_ph(state.site, midpoint)).background
        if color != current_color:
            segments.append(
                PkaPlotSegment(
                    start_ph=current_start,
                    end_ph=start_ph,
                    color=current_color,
                )
            )
            current_start = start_ph
            current_color = color
        ph = end_ph

    segments.append(
        PkaPlotSegment(
            start_ph=current_start,
            end_ph=PH_MAX,
            color=current_color,
        )
    )
    return tuple(segments)


def _clip_to_plot(value: float, inset: float = 0.0) -> float:
    lower = PH_MIN + inset
    upper = PH_MAX - inset
    if value < lower:
        return lower
    if value > upper:
        return upper
    return value
