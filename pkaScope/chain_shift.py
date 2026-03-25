from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import ceil
from typing import Iterable

from matplotlib import pyplot as plt
from matplotlib.figure import Figure

from .pdb_utils import ResidueKey, atom_residue_key, load_pdb_universe, write_pdb_atoms
from .propka_analysis import TitrationSite

SiteMatchKey = tuple[int, str, str]


@dataclass(frozen=True)
class ChainPkaShift:
    label: str
    residue_type: str
    chain_id: str
    residue_number: int
    insertion_code: str
    complex_pka: float
    monomer_pka: float
    delta_pka: float
    absolute_delta_pka: float

    @property
    def residue_key(self) -> ResidueKey:
        return (self.chain_id, self.residue_number, self.insertion_code)

    @property
    def residue_label(self) -> str:
        insertion_code = self.insertion_code.strip()
        residue_token = (
            f"{self.residue_number}{insertion_code}"
            if insertion_code
            else str(self.residue_number)
        )
        return f"{self.chain_id}:{residue_token}-{self.residue_type}"


@dataclass(frozen=True)
class ChainPkaComparison:
    chain_id: str
    shifts: tuple[ChainPkaShift, ...]
    unmatched_complex_labels: tuple[str, ...]
    unmatched_monomer_labels: tuple[str, ...]

    @property
    def max_absolute_delta(self) -> float:
        if not self.shifts:
            return 0.0
        return max(shift.absolute_delta_pka for shift in self.shifts)


def list_chain_ids(pdb_text: str) -> tuple[str, ...]:
    universe = load_pdb_universe(pdb_text)
    ordered_chain_ids: dict[str, None] = {}
    for atom in universe.atoms:
        chain_id = atom_residue_key(atom)[0]
        ordered_chain_ids.setdefault(chain_id, None)
    return tuple(ordered_chain_ids.keys())


def extract_chain_pdb(pdb_text: str, chain_id: str) -> str:
    universe = load_pdb_universe(pdb_text)
    normalized_chain_id = _normalize_chain_id(chain_id)
    selected_atom_indices = [
        atom.index
        for atom in universe.atoms
        if atom_residue_key(atom)[0] == normalized_chain_id
    ]
    if not selected_atom_indices:
        raise ValueError(f"Chain '{normalized_chain_id}' has no coordinate records.")
    return write_pdb_atoms(universe.atoms[selected_atom_indices])


def compare_chain_pkas(
    *,
    chain_id: str,
    complex_sites: Iterable[TitrationSite],
    monomer_sites: Iterable[TitrationSite],
) -> ChainPkaComparison:
    normalized_chain_id = _normalize_chain_id(chain_id)
    complex_chain_sites = tuple(
        site for site in complex_sites if _normalize_chain_id(site.chain_id) == normalized_chain_id
    )
    monomer_chain_sites = tuple(
        site for site in monomer_sites if _normalize_chain_id(site.chain_id) == normalized_chain_id
    )
    if not monomer_chain_sites:
        monomer_chain_sites = tuple(monomer_sites)

    remaining_monomer_by_site_key = _group_sites_by_site_key(monomer_chain_sites)
    remaining_monomer_by_match_key = _group_sites_by_match_key(monomer_chain_sites)

    shifts: list[ChainPkaShift] = []
    unmatched_complex_labels: list[str] = []

    for complex_site in sorted(complex_chain_sites, key=_site_sort_key):
        monomer_site = _pop_site(remaining_monomer_by_site_key, complex_site.site_key)
        if monomer_site is None:
            monomer_site = _pop_site(
                remaining_monomer_by_match_key,
                _site_match_key(complex_site),
            )
            if monomer_site is not None:
                _discard_site(
                    remaining_monomer_by_site_key,
                    monomer_site.site_key,
                    monomer_site,
                )
        if monomer_site is None:
            unmatched_complex_labels.append(complex_site.label)
            continue
        delta_pka = complex_site.pka - monomer_site.pka
        shifts.append(
            ChainPkaShift(
                label=complex_site.label,
                residue_type=complex_site.residue_type,
                chain_id=complex_site.chain_id,
                residue_number=complex_site.residue_number,
                insertion_code=complex_site.insertion_code,
                complex_pka=complex_site.pka,
                monomer_pka=monomer_site.pka,
                delta_pka=delta_pka,
                absolute_delta_pka=abs(delta_pka),
            )
        )

    unmatched_monomer_labels = sorted(
        site.label
        for remaining_sites in remaining_monomer_by_site_key.values()
        for site in remaining_sites
    )

    return ChainPkaComparison(
        chain_id=normalized_chain_id,
        shifts=tuple(sorted(shifts, key=_shift_sort_key)),
        unmatched_complex_labels=tuple(sorted(unmatched_complex_labels)),
        unmatched_monomer_labels=tuple(unmatched_monomer_labels),
    )


def top_shift_residue_keys(
    shifts: Iterable[ChainPkaShift],
    *,
    minimum_absolute_delta: float = 0.5,
    limit: int = 8,
) -> tuple[ResidueKey, ...]:
    ranked_shifts = sorted(
        (shift for shift in shifts if shift.absolute_delta_pka >= minimum_absolute_delta),
        key=lambda shift: (
            -shift.absolute_delta_pka,
            shift.chain_id,
            shift.residue_number,
            shift.insertion_code,
            shift.label,
        ),
    )
    return tuple(shift.residue_key for shift in ranked_shifts[:limit])


def create_chain_shift_plot_figure(
    shifts: Iterable[ChainPkaShift],
    current_ph: float | None = None,
) -> Figure:
    ordered_shifts = list(shifts)
    if not ordered_shifts:
        figure, axis = plt.subplots(figsize=(10.0, 3.0))
        axis.text(
            0.5,
            0.5,
            "No matched titratable sites were available for this chain.",
            ha="center",
            va="center",
            fontsize=12,
        )
        axis.axis("off")
        figure.tight_layout()
        return figure

    complex_values = [shift.complex_pka for shift in ordered_shifts]
    monomer_values = [shift.monomer_pka for shift in ordered_shifts]
    delta_values = [shift.delta_pka for shift in ordered_shifts]
    labels = [shift.residue_label for shift in ordered_shifts]
    x_positions = list(range(len(ordered_shifts)))

    figure_width = max(10.0, (len(ordered_shifts) * 0.4) + 3.5)
    figure, (pka_axis, delta_axis) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(figure_width, 3.6),
        gridspec_kw={"height_ratios": (3.2, 1.5)},
    )

    pka_axis.plot(
        x_positions,
        complex_values,
        marker="o",
        markersize=4.4,
        linewidth=1.8,
        color="#1d4ed8",
        label="Complex",
    )
    pka_axis.plot(
        x_positions,
        monomer_values,
        marker="s",
        markersize=4.4,
        linewidth=1.8,
        color="#ea580c",
        label="Monomer",
    )
    if current_ph is not None:
        pka_axis.axhline(
            current_ph,
            color="#111827",
            linewidth=1.2,
            linestyle="--",
            alpha=0.85,
            label=f"Solution pH {current_ph:.1f}",
        )
    for x_position, complex_pka, monomer_pka in zip(
        x_positions, complex_values, monomer_values
    ):
        pka_axis.vlines(
            x_position,
            ymin=min(complex_pka, monomer_pka),
            ymax=max(complex_pka, monomer_pka),
            color="#6b7280",
            linewidth=1.2,
            alpha=0.55,
            zorder=1,
        )
    pka_axis.set_ylabel("pKa")
    pka_axis.grid(axis="y", color="#9ca3af", alpha=0.35, linewidth=0.8)
    pka_axis.legend(loc="upper left", frameon=False)
    pka_axis.spines["top"].set_visible(False)
    pka_axis.spines["right"].set_visible(False)

    delta_axis.bar(
        x_positions,
        delta_values,
        color=[_delta_color(value) for value in delta_values],
        alpha=0.85,
        width=0.72,
    )
    delta_axis.axhline(0.0, color="#111827", linewidth=1.1, alpha=0.85)
    delta_axis.grid(axis="y", color="#9ca3af", alpha=0.35, linewidth=0.8)
    delta_axis.set_ylabel("Delta pKa")
    delta_axis.set_xlabel("Chain residues")
    delta_axis.spines["top"].set_visible(False)
    delta_axis.spines["right"].set_visible(False)

    tick_step = max(1, ceil(len(labels) / 30))
    tick_positions = x_positions[::tick_step]
    tick_labels = labels[::tick_step]
    delta_axis.set_xticks(tick_positions)
    delta_axis.set_xticklabels(tick_labels, rotation=70, ha="right", fontsize=8)

    figure.tight_layout()
    return figure


def _normalize_chain_id(chain_id: str) -> str:
    return chain_id.strip() or "?"


def _group_sites_by_site_key(
    sites: Iterable[TitrationSite],
) -> dict[tuple[str, int, str, str], list[TitrationSite]]:
    grouped_sites: dict[tuple[str, int, str, str], list[TitrationSite]] = defaultdict(list)
    for site in sites:
        grouped_sites[site.site_key].append(site)
    return grouped_sites


def _group_sites_by_match_key(
    sites: Iterable[TitrationSite],
) -> dict[SiteMatchKey, list[TitrationSite]]:
    grouped_sites: dict[SiteMatchKey, list[TitrationSite]] = defaultdict(list)
    for site in sites:
        grouped_sites[_site_match_key(site)].append(site)
    return grouped_sites


def _site_match_key(site: TitrationSite) -> SiteMatchKey:
    return (site.residue_number, site.insertion_code, site.residue_type)


def _pop_site[T](grouped_sites: dict[T, list[TitrationSite]], key: T) -> TitrationSite | None:
    grouped = grouped_sites.get(key)
    if not grouped:
        return None
    site = grouped.pop(0)
    if not grouped:
        del grouped_sites[key]
    return site


def _discard_site[T](
    grouped_sites: dict[T, list[TitrationSite]],
    key: T,
    site: TitrationSite,
) -> None:
    grouped = grouped_sites.get(key)
    if not grouped:
        return
    for index, candidate_site in enumerate(grouped):
        if candidate_site == site:
            del grouped[index]
            if not grouped:
                del grouped_sites[key]
            return


def _site_sort_key(site: TitrationSite) -> tuple[int, str, str]:
    return (site.residue_number, site.insertion_code, site.label)


def _shift_sort_key(shift: ChainPkaShift) -> tuple[int, str, str]:
    return (shift.residue_number, shift.insertion_code, shift.label)


def _delta_color(delta_pka: float) -> str:
    if delta_pka >= 0.5:
        return "#16a34a"
    if delta_pka <= -0.5:
        return "#dc2626"
    return "#6b7280"
