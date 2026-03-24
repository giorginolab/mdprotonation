from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ChargeColorBand:
    background: str
    foreground: str


@dataclass(frozen=True)
class ChargeLegendBand:
    label: str
    color: str


def charge_color_band(charge: float) -> ChargeColorBand:
    if charge <= -0.75:
        return ChargeColorBand(background="#c62828", foreground="#ffffff")
    if charge <= -0.25:
        return ChargeColorBand(background="#f6c1d6", foreground="#3f1025")
    if charge < 0.25:
        return ChargeColorBand(background="#ffffff", foreground="#111827")
    if charge < 0.75:
        return ChargeColorBand(background="#6f82be", foreground="#ffffff")
    return ChargeColorBand(background="#0d47a1", foreground="#ffffff")


def charge_legend_bands() -> tuple[ChargeLegendBand, ...]:
    return (
        ChargeLegendBand("<= -0.75", charge_color_band(-1.0).background),
        ChargeLegendBand("(-0.75, -0.25]", charge_color_band(-0.5).background),
        ChargeLegendBand("(-0.25, +0.25)", charge_color_band(0.0).background),
        ChargeLegendBand("[+0.25, +0.75)", charge_color_band(0.5).background),
        ChargeLegendBand(">= +0.75", charge_color_band(1.0).background),
    )
