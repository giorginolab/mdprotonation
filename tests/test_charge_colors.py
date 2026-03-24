from __future__ import annotations

from pkaScope.charge_colors import charge_color_band


def test_charge_color_band_uses_protonation_threshold_boundaries() -> None:
    # Acidic side: <10% prot (charge < -0.9) is dark red.
    assert charge_color_band(-0.91).background == "#c62828"
    # Acidic transition window includes the -0.9 boundary.
    assert charge_color_band(-0.90).background == "#f6c1d6"
    assert charge_color_band(-0.11).background == "#f6c1d6"
    # Acidic >90% prot / basic <10% prot is white.
    assert charge_color_band(-0.10).background == "#ffffff"
    assert charge_color_band(0.09).background == "#ffffff"
    # Basic transition window includes the +0.1 boundary.
    assert charge_color_band(0.10).background == "#6f82be"
    assert charge_color_band(0.89).background == "#6f82be"
    # Basic >90% prot starts at +0.9.
    assert charge_color_band(0.90).background == "#0d47a1"
