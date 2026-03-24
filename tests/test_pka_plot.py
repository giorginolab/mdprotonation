from __future__ import annotations

from mdprotonation.pka_plot import (
    build_pka_plot_rows,
    create_pka_plot_figure,
    create_pka_comparison_plot_figure,
)
from mdprotonation.charge_colors import charge_color_band
from mdprotonation.propka_analysis import TitrationSite
from mdprotonation.protonation import SiteState


def make_site_state(
    *,
    residue_type: str,
    chain_id: str,
    residue_number: int,
    insertion_code: str = "",
    pka: float,
    model_pka: float | None = None,
    charged_state_charge: float,
    dominant_state: str = "Neutral",
) -> SiteState:
    site = TitrationSite(
        label=f"{residue_type} {residue_number:>3} {chain_id}",
        residue_type=residue_type,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code=insertion_code,
        pka=pka,
        model_pka=model_pka if model_pka is not None else pka,
        charged_state_charge=charged_state_charge,
        buried_fraction=0.0,
        interactions=(),
    )
    return SiteState(
        site=site,
        ph=7.0,
        protonated_fraction=0.5,
        deprotonated_fraction=0.5,
        current_charge=0.0,
        transition_score=1.0,
        dominant_state=dominant_state,
    )


def test_build_pka_plot_rows_formats_labels_and_charge_ranges() -> None:
    positive_state = make_site_state(
        residue_type="LYS",
        chain_id="B",
        residue_number=44,
        insertion_code="A",
        pka=10.5,
        charged_state_charge=1.0,
    )
    negative_state = make_site_state(
        residue_type="ASP",
        chain_id="A",
        residue_number=12,
        pka=3.7,
        charged_state_charge=-1.0,
        dominant_state="Transitioning",
    )

    rows = build_pka_plot_rows([positive_state, negative_state])

    assert [row.label for row in rows] == ["⚠️ A:12-ASP", "B:44A-LYS"]
    assert rows[0].charge_segments[0].color == charge_color_band(0.0).background
    assert rows[0].charge_segments[-1].color == charge_color_band(-1.0).background
    assert rows[1].charge_segments[0].color == charge_color_band(1.0).background
    assert rows[1].charge_segments[-1].color == charge_color_band(0.0).background
    assert any(
        segment.color == charge_color_band(-0.5).background
        for segment in rows[0].charge_segments
    )
    assert any(
        segment.color == charge_color_band(0.5).background
        for segment in rows[1].charge_segments
    )


def test_build_pka_plot_rows_clips_out_of_range_markers() -> None:
    highly_basic_state = make_site_state(
        residue_type="ARG",
        chain_id="A",
        residue_number=101,
        pka=18.2,
        charged_state_charge=1.0,
    )
    strongly_acidic_state = make_site_state(
        residue_type="ASP",
        chain_id="A",
        residue_number=18,
        pka=-1.6,
        charged_state_charge=-1.0,
    )

    rows = build_pka_plot_rows([highly_basic_state, strongly_acidic_state])

    assert [row.label for row in rows] == ["A:18-ASP", "A:101-ARG"]
    assert rows[0].pka_marker == 0.22
    assert rows[0].charge_segments[0].start_ph == 0.0
    assert rows[0].charge_segments[-1].end_ph == 14.0
    assert rows[1].pka_marker == 13.78
    assert rows[1].charge_segments[0].start_ph == 0.0
    assert rows[1].charge_segments[-1].end_ph == 14.0


def test_create_pka_comparison_plot_figure_runs_without_error() -> None:
    states = [
        make_site_state(
            residue_type="LYS",
            chain_id="A",
            residue_number=1,
            pka=10.5,
            model_pka=10.0,
            charged_state_charge=1.0,
        ),
        make_site_state(
            residue_type="ASP",
            chain_id="A",
            residue_number=2,
            pka=3.7,
            model_pka=4.0,
            charged_state_charge=-1.0,
        ),
    ]

    figure = create_pka_comparison_plot_figure(states, current_ph=7.4)
    assert figure is not None
    assert len(figure.axes) == 2


def test_create_pka_plot_figure_returns_plotly_figure() -> None:
    states = [
        make_site_state(
            residue_type="LYS",
            chain_id="A",
            residue_number=1,
            pka=10.5,
            model_pka=10.0,
            charged_state_charge=1.0,
        ),
        make_site_state(
            residue_type="ASP",
            chain_id="A",
            residue_number=2,
            pka=3.7,
            model_pka=4.0,
            charged_state_charge=-1.0,
            dominant_state="Transitioning",
        ),
    ]

    figure = create_pka_plot_figure(states, current_ph=7.4)

    assert figure is not None
    assert len(figure.data) >= 3
    assert figure.layout.xaxis.title.text == "pKa"
    assert figure.layout.yaxis.autorange == "reversed"
