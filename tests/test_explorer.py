from __future__ import annotations

from mdprotonation.panes.explorer import _site_column_label
from mdprotonation.panes.explorer import (
    SITE_TABLE_MAX_HEIGHT,
    SITE_TABLE_MIN_HEIGHT,
    _site_table_height,
)
from mdprotonation.propka_analysis import TitrationSite
from mdprotonation.protonation import SiteState


def test_site_table_height_enforces_minimum_for_small_tables() -> None:
    assert _site_table_height(0) == SITE_TABLE_MIN_HEIGHT
    assert _site_table_height(3) == SITE_TABLE_MIN_HEIGHT


def test_site_table_height_caps_large_tables_at_viewer_height() -> None:
    assert _site_table_height(100) == SITE_TABLE_MAX_HEIGHT


def test_site_column_label_uses_chain_and_natural_residue_token() -> None:
    site = TitrationSite(
        label="ASP  12 A",
        residue_type="ASP",
        chain_id="A",
        residue_number=12,
        insertion_code="B",
        pka=4.0,
        model_pka=4.0,
        charged_state_charge=-1.0,
        buried_fraction=0.1,
        interactions=(),
    )
    state = SiteState(
        site=site,
        ph=7.0,
        protonated_fraction=0.2,
        deprotonated_fraction=0.8,
        current_charge=-0.8,
        transition_score=0.4,
        dominant_state="Anionic",
    )

    assert _site_column_label(state) == "A:   12B-ASP"
