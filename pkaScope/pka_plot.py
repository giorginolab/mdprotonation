from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import ceil

import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
PLOT_FONT_FAMILY = "Avenir Next, Helvetica Neue, Arial, sans-serif"


@dataclass(frozen=True)
class PkaPlotRow:
    label: str
    site_label: str
    pka: float
    pka_marker: float
    charge_segments: tuple["PkaPlotSegment", ...]
    is_transitioning: bool
    dominant_state: str


@dataclass(frozen=True)
class PkaPlotSegment:
    start_ph: float
    end_ph: float
    color: str
    representative_charge: float


def build_pka_plot_rows(site_states: list[SiteState]) -> list[PkaPlotRow]:
    rows: list[PkaPlotRow] = []
    sorted_states = sorted(site_states, key=site_state_sort_key)
    for state in sorted_states:
        rows.append(
            PkaPlotRow(
                label=_format_plot_label(state),
                site_label=state.site.label,
                pka=state.site.pka,
                pka_marker=_clip_to_plot(state.site.pka, inset=0.22),
                charge_segments=_build_charge_segments(state),
                is_transitioning=state.dominant_state == "Transitioning",
                dominant_state=state.dominant_state,
            )
        )
    return rows


def create_pka_plot_figure(site_states: list[SiteState], current_ph: float) -> go.Figure:
    rows = build_pka_plot_rows(site_states)
    figure = go.Figure()
    segment_buckets: dict[str, list[tuple[PkaPlotRow, PkaPlotSegment]]] = defaultdict(list)
    category_order = [row.label for row in rows]

    for row in rows:
        for segment in row.charge_segments:
            segment_buckets[segment.color].append((row, segment))

    for band in charge_legend_bands():
        entries = segment_buckets.get(band.color, [])
        if not entries:
            continue
        figure.add_trace(
            go.Bar(
                name=band.label,
                orientation="h",
                x=[segment.end_ph - segment.start_ph for _, segment in entries],
                base=[segment.start_ph for _, segment in entries],
                y=[row.label for row, _ in entries],
                width=0.82,
                marker={
                    "color": band.color,
                    "line": {
                        "color": LANE_EDGE_COLOR if band.color == LANE_COLOR else band.color,
                        "width": 0.7 if band.color == LANE_COLOR else 0.0,
                    },
                },
                customdata=[
                    [
                        row.site_label,
                        row.pka,
                        segment.start_ph,
                        segment.end_ph,
                        segment.representative_charge,
                        row.dominant_state,
                    ]
                    for row, segment in entries
                ],
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "pKa: %{customdata[1]:.2f}<br>"
                    "pH range: %{customdata[2]:.2f} to %{customdata[3]:.2f}<br>"
                    "Representative charge: %{customdata[4]:.2f}<br>"
                    "State at selected pH: %{customdata[5]}<extra></extra>"
                ),
            )
        )

    for row in rows:
        figure.add_annotation(
            x=row.pka_marker,
            y=row.label,
            text=f"{row.pka:.2f}",
            showarrow=False,
            font={"size": 11, "color": "#1f2937", "family": PLOT_FONT_FAMILY},
            bgcolor="rgba(248, 250, 252, 0.9)",
            opacity=0.95,
            borderpad=2,
        )

    for boundary in (current_ph - 1.0, current_ph + 1.0):
        if PH_MIN <= boundary <= PH_MAX:
            figure.add_vline(
                x=boundary,
                line_color=TRANSITION_WINDOW_COLOR,
                line_width=1.6,
                line_dash="dash",
                opacity=0.55,
            )
    figure.add_vline(
        x=current_ph,
        line_color=CURRENT_PH_COLOR,
        line_width=2.4,
        opacity=0.92,
    )

    figure.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"Selected pH {current_ph:.1f}",
            line={"color": CURRENT_PH_COLOR, "width": 2.4},
            hoverinfo="skip",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=TRANSITION_WINDOW_LABEL,
            line={"color": TRANSITION_WINDOW_COLOR, "width": 1.6, "dash": "dash"},
            hoverinfo="skip",
        )
    )

    figure.update_layout(
        barmode="overlay",
        bargap=0.18,
        height=max(540, 170 + (len(rows) * 28)),
        margin={"l": 24, "r": 24, "t": 132, "b": 28},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#fcfbf8",
        font={"family": PLOT_FONT_FAMILY, "size": 14, "color": "#111827"},
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.1,
            "xanchor": "left",
            "x": 0.0,
            "font": {"size": 13},
        },
        hovermode="closest",
        hoverlabel={
            "bgcolor": "#ffffff",
            "bordercolor": "#d1d5db",
            "font": {"family": PLOT_FONT_FAMILY, "size": 13, "color": "#111827"},
        },
    )
    figure.update_xaxes(
        range=[PH_MIN, PH_MAX],
        side="top",
        title_text="pKa",
        title_font={"size": 20},
        dtick=2,
        tick0=0,
        tickfont={"size": 14},
        showgrid=True,
        gridcolor="rgba(55, 65, 81, 0.24)",
        gridwidth=1.1,
        showline=True,
        linecolor="#9ca3af",
        linewidth=1.0,
        zeroline=False,
        minor={"dtick": 1, "showgrid": True, "gridcolor": "rgba(107, 114, 128, 0.18)"},
    )
    figure.update_yaxes(
        autorange="reversed",
        categoryorder="array",
        categoryarray=category_order,
        showgrid=False,
        ticks="",
        tickfont={"size": 14},
        automargin=True,
    )
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
    current_charge = site_charge_at_ph(state.site, PH_MIN + (PH_SEGMENT_STEP / 2.0))
    current_color = charge_color_band(current_charge).background

    ph = PH_MIN
    while ph < PH_MAX:
        start_ph = ph
        end_ph = min(ph + PH_SEGMENT_STEP, PH_MAX)
        midpoint = start_ph + ((end_ph - start_ph) / 2.0)
        midpoint_charge = site_charge_at_ph(state.site, midpoint)
        color = charge_color_band(midpoint_charge).background
        if color != current_color:
            segments.append(
                PkaPlotSegment(
                    start_ph=current_start,
                    end_ph=start_ph,
                    color=current_color,
                    representative_charge=current_charge,
                )
            )
            current_start = start_ph
            current_color = color
            current_charge = midpoint_charge
        ph = end_ph

    segments.append(
        PkaPlotSegment(
            start_ph=current_start,
            end_ph=PH_MAX,
            color=current_color,
            representative_charge=current_charge,
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


def create_pka_comparison_plot_figure(site_states: list[SiteState], current_ph: float) -> go.Figure:
    sorted_states = sorted(site_states, key=site_state_sort_key)
    if not sorted_states:
        figure = go.Figure()
        figure.add_annotation(
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            text="No titratable sites were found.",
            showarrow=False,
            font={"size": 12, "family": PLOT_FONT_FAMILY, "color": "#111827"},
        )
        figure.update_layout(
            height=320,
            margin={"l": 24, "r": 24, "t": 24, "b": 24},
            xaxis={"visible": False},
            yaxis={"visible": False},
            paper_bgcolor="#ffffff",
            plot_bgcolor="#ffffff",
        )
        return figure

    computed_values = [state.site.pka for state in sorted_states]
    model_values = [state.site.model_pka for state in sorted_states]
    delta_values = [pka - model for pka, model in zip(computed_values, model_values)]
    labels = [_format_plot_label(state) for state in sorted_states]
    x_positions = list(range(len(sorted_states)))

    figure = make_subplots(
        2,
        1,
        shared_xaxes=True,
        row_heights=[0.68, 0.32],
        vertical_spacing=0.06,
    )

    figure.add_trace(
        go.Scatter(
            x=x_positions,
            y=computed_values,
            mode="lines+markers",
            marker={"symbol": "circle", "size": 8},
            line={"width": 2.5, "color": "#1d4ed8"},
            name="Computed pKa",
            customdata=[[labels[index]] for index in x_positions],
            hovertemplate="<b>%{customdata[0]}</b><br>Computed pKa: %{y:.2f}<extra></extra>",
        ),
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=x_positions,
            y=model_values,
            mode="lines+markers",
            marker={"symbol": "square", "size": 8},
            line={"width": 2.5, "color": "#ea580c"},
            name="Model pKa",
            customdata=[[labels[index]] for index in x_positions],
            hovertemplate="<b>%{customdata[0]}</b><br>Model pKa: %{y:.2f}<extra></extra>",
        ),
        row=1,
        col=1,
    )
    figure.add_hline(
        y=current_ph,
        line_color="#111827",
        line_width=1.8,
        line_dash="dash",
        opacity=0.85,
        row=1,
        col=1,
    )
    figure.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode="lines",
            name=f"Solution pH {current_ph:.1f}",
            line={"color": "#111827", "width": 1.8, "dash": "dash"},
            hoverinfo="skip",
        ),
        row=1,
        col=1,
    )
    for x_pos, comp_pka, mod_pka in zip(x_positions, computed_values, model_values):
        figure.add_trace(
            go.Scatter(
                x=[x_pos, x_pos],
                y=[min(comp_pka, mod_pka), max(comp_pka, mod_pka)],
                mode="lines",
                showlegend=False,
                line={"color": "#6b7280", "width": 1.8},
                opacity=0.55,
                hoverinfo="skip",
            ),
            row=1,
            col=1,
        )

    figure.add_trace(
        go.Bar(
            x=x_positions,
            y=delta_values,
            marker={"color": [_delta_color(value) for value in delta_values]},
            opacity=0.85,
            width=0.72,
            showlegend=False,
            customdata=[[labels[index]] for index in x_positions],
            hovertemplate="<b>%{customdata[0]}</b><br>Delta pKa: %{y:+.2f}<extra></extra>",
        ),
        row=2,
        col=1,
    )
    figure.add_hline(
        y=0.0,
        line_color="#111827",
        line_width=1.5,
        opacity=0.85,
        row=2,
        col=1,
    )

    tick_step = max(1, ceil(len(labels) / 30))
    tick_positions = x_positions[::tick_step]
    tick_labels = labels[::tick_step]
    figure.update_xaxes(
        tickmode="array",
        tickvals=tick_positions,
        ticktext=tick_labels,
        tickangle=70,
        title_text="Residues",
        tickfont={"size": 12},
        row=2,
        col=1,
    )
    figure.update_yaxes(
        title_text="pKa",
        tickfont={"size": 14},
        showgrid=True,
        gridcolor="rgba(156, 163, 175, 0.35)",
        gridwidth=0.8,
        row=1,
        col=1,
    )
    figure.update_yaxes(
        title_text="Delta pKa",
        tickfont={"size": 14},
        showgrid=True,
        gridcolor="rgba(156, 163, 175, 0.35)",
        gridwidth=0.8,
        row=2,
        col=1,
    )
    figure.update_layout(
        height=max(560, 220 + (len(sorted_states) * 16)),
        margin={"l": 36, "r": 24, "t": 44, "b": 24},
        font={"family": PLOT_FONT_FAMILY, "size": 14, "color": "#111827"},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        legend={
            "orientation": "h",
            "yanchor": "bottom",
            "y": 1.02,
            "xanchor": "left",
            "x": 0.0,
        },
        hovermode="x unified",
    )

    return figure


def _delta_color(delta_pka: float) -> str:
    if delta_pka >= 0.5:
        return "#16a34a"
    if delta_pka <= -0.5:
        return "#dc2626"
    return "#6b7280"
