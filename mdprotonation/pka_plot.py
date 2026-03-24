from __future__ import annotations

from dataclasses import dataclass

from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .protonation import SiteState

PH_MIN = 0.0
PH_MAX = 14.0
LANE_COLOR = "#b8b8b8"
POSITIVE_RANGE_COLOR = "#ff6b4a"
NEGATIVE_RANGE_COLOR = "#74b3d8"
CURRENT_PH_COLOR = "#3f3f46"


@dataclass(frozen=True)
class PkaPlotRow:
    label: str
    pka: float
    pka_marker: float
    charged_range_start: float
    charged_range_end: float
    charged_range_color: str
    is_transitioning: bool


def build_pka_plot_rows(site_states: list[SiteState]) -> list[PkaPlotRow]:
    rows: list[PkaPlotRow] = []
    sorted_states = sorted(
        site_states,
        key=lambda state: (
            state.site.chain_id,
            state.site.residue_number,
            state.site.insertion_code,
            state.site.label,
        ),
    )
    for state in sorted_states:
        site = state.site
        if site.charged_state_charge > 0:
            charged_range_start = PH_MIN
            charged_range_end = _clip_to_plot(site.pka)
            charged_range_color = POSITIVE_RANGE_COLOR
        else:
            charged_range_start = _clip_to_plot(site.pka)
            charged_range_end = PH_MAX
            charged_range_color = NEGATIVE_RANGE_COLOR

        rows.append(
            PkaPlotRow(
                label=_format_plot_label(state),
                pka=site.pka,
                pka_marker=_clip_to_plot(site.pka, inset=0.22),
                charged_range_start=charged_range_start,
                charged_range_end=charged_range_end,
                charged_range_color=charged_range_color,
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
            color=LANE_COLOR,
            edgecolor="white",
            linewidth=0.35,
            zorder=1,
        )
        charged_width = row.charged_range_end - row.charged_range_start
        if charged_width > 0:
            axis.barh(
                index,
                charged_width,
                left=row.charged_range_start,
                height=0.82,
                color=row.charged_range_color,
                edgecolor="white",
                linewidth=0.35,
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
        Patch(facecolor=POSITIVE_RANGE_COLOR, label="Positive charged range"),
        Patch(facecolor=NEGATIVE_RANGE_COLOR, label="Negative charged range"),
        Patch(facecolor=LANE_COLOR, label="Outside charged range"),
        Line2D(
            [0],
            [0],
            color=CURRENT_PH_COLOR,
            linewidth=2.4,
            label=f"Selected pH {current_ph:.1f}",
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
    prefix = "(!) " if state.dominant_state == "Transitioning" else ""
    return f"{prefix}{site.chain_id}:{residue_token}-{site.residue_type}"


def _clip_to_plot(value: float, inset: float = 0.0) -> float:
    lower = PH_MIN + inset
    upper = PH_MAX - inset
    if value < lower:
        return lower
    if value > upper:
        return upper
    return value
