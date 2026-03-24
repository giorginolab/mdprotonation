from __future__ import annotations

from dataclasses import dataclass

LOW_PROTONATION_THRESHOLD = 0.1
HIGH_PROTONATION_THRESHOLD = 0.9

# Charge-space equivalents for the protonation thresholds:
# acidic charge = -(1 - prot), basic charge = prot
LOW_CHARGE_THRESHOLD = -HIGH_PROTONATION_THRESHOLD
HIGH_CHARGE_THRESHOLD = -LOW_PROTONATION_THRESHOLD
NEUTRAL_BASIC_THRESHOLD = LOW_PROTONATION_THRESHOLD
TRANSITION_BASIC_THRESHOLD = HIGH_PROTONATION_THRESHOLD


@dataclass(frozen=True)
class ChargeColorBand:
    background: str
    foreground: str


@dataclass(frozen=True)
class ChargeLegendBand:
    label: str
    color: str


def charge_color_band(charge: float) -> ChargeColorBand:
    if charge < LOW_CHARGE_THRESHOLD:
        return ChargeColorBand(background="#c62828", foreground="#ffffff")
    if charge < HIGH_CHARGE_THRESHOLD:
        return ChargeColorBand(background="#f6c1d6", foreground="#3f1025")
    if charge < NEUTRAL_BASIC_THRESHOLD:
        return ChargeColorBand(background="#ffffff", foreground="#111827")
    if charge < TRANSITION_BASIC_THRESHOLD:
        return ChargeColorBand(background="#6f82be", foreground="#ffffff")
    return ChargeColorBand(background="#0d47a1", foreground="#ffffff")


def charge_legend_bands() -> tuple[ChargeLegendBand, ...]:
    return (
        ChargeLegendBand("Acidic: <10% prot", charge_color_band(-1.0).background),
        ChargeLegendBand("Acidic: 10-90% prot", charge_color_band(-0.5).background),
        ChargeLegendBand(
            "Acidic >90% / Basic <10% prot",
            charge_color_band(0.0).background,
        ),
        ChargeLegendBand("Basic: 10-90% prot", charge_color_band(0.5).background),
        ChargeLegendBand("Basic: >90% prot", charge_color_band(1.0).background),
    )
