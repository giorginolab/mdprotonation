from __future__ import annotations

import pandas as pd

from .charge_colors import charge_color_band

TRANSITION_MARKER = "⚠️"
LANE_COLOR = charge_color_band(0.0).background
NEGATIVE_CHARGE_HIGH_COLOR = charge_color_band(-1.0).background
NEGATIVE_CHARGE_MID_COLOR = charge_color_band(-0.5).background
NEUTRAL_CHARGE_COLOR = charge_color_band(0.0).background
POSITIVE_CHARGE_MID_COLOR = charge_color_band(0.5).background
POSITIVE_CHARGE_HIGH_COLOR = charge_color_band(1.0).background
CURRENT_PH_COLOR = "#3f3f46"
TRANSITION_WINDOW_COLOR = "#6b7280"
NEGATIVE_RANGE_COLOR = NEGATIVE_CHARGE_HIGH_COLOR
POSITIVE_RANGE_COLOR = POSITIVE_CHARGE_HIGH_COLOR

CHARGE_SCALE_CAPTION = (
    "Row colors by protonation class: acidic <10% prot red, acidic 10-90% prot pink, "
    "acidic >90% or basic <10% prot white, basic 10-90% prot slate blue, "
    "basic >90% prot deep blue."
)
PH7_STATE_CAPTION = f"`{TRANSITION_MARKER}` marks sites not in their pH 7 dominant state."
VIEWER_ENCODING_CAPTION = (
    "Viewer encoding: occupancy stores average protonated fraction per "
    "titratable residue and B-factor stores transition intensity x100."
)
PKA_PLOT_CAPTION = (
    "The pKa lanes use the same five acidic/basic protonation bins as the Explorer "
    "table: red, pink, white, slate blue, and deep blue. The solid vertical line "
    "marks the selected pH, dashed lines mark `pH ± 1`, and "
    f"`{TRANSITION_MARKER}` marks sites transitioning at the selected pH."
)
POSITIVE_RANGE_LABEL = "Positive charged range (electron-deficient)"
NEGATIVE_RANGE_LABEL = "Negative charged range (electron-rich)"
TRANSITION_WINDOW_LABEL = "pH ± 1 transition window"


def charge_palette(charge: float) -> tuple[str, str]:
    color_band = charge_color_band(charge)
    return (color_band.background, color_band.foreground)


def style_charge_row(row: pd.Series) -> list[str]:
    background, foreground = charge_palette(float(row["Charge"]))
    return [f"background-color: {background}; color: {foreground}"] * len(row)
