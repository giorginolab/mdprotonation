from __future__ import annotations

from mdprotonation.panes.explorer import (
    SITE_TABLE_MAX_HEIGHT,
    SITE_TABLE_MIN_HEIGHT,
    _site_table_height,
)


def test_site_table_height_enforces_minimum_for_small_tables() -> None:
    assert _site_table_height(0) == SITE_TABLE_MIN_HEIGHT
    assert _site_table_height(3) == SITE_TABLE_MIN_HEIGHT


def test_site_table_height_caps_large_tables_at_viewer_height() -> None:
    assert _site_table_height(100) == SITE_TABLE_MAX_HEIGHT
