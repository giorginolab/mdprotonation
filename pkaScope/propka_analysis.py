from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Iterable

import propka.output
import propka.run


@dataclass(frozen=True)
class SiteInteraction:
    partner_label: str
    partner_residue_type: str
    chain_id: str
    residue_number: int
    insertion_code: str
    interaction_kind: str
    contribution: float


@dataclass(frozen=True)
class TitrationSite:
    label: str
    residue_type: str
    chain_id: str
    residue_number: int
    insertion_code: str
    pka: float
    model_pka: float
    charged_state_charge: float
    buried_fraction: float
    interactions: tuple[SiteInteraction, ...]

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain_id, self.residue_number, self.insertion_code)

    @property
    def site_key(self) -> tuple[str, int, str, str]:
        return (*self.residue_key, self.label)


@dataclass(frozen=True)
class ChargePoint:
    ph: float
    unfolded_charge: float
    folded_charge: float


@dataclass(frozen=True)
class FoldingEnergyPoint:
    ph: float
    folding_free_energy: float


@dataclass(frozen=True)
class PropkaAnalysis:
    source_name: str
    structure_name: str
    pdb_text: str
    titration_sites: tuple[TitrationSite, ...]
    charge_profile: tuple[ChargePoint, ...]
    folding_profile: tuple[FoldingEnergyPoint, ...]
    folding_optimum_ph: float | None
    folding_optimum_energy: float | None
    summary_text: str
    determinants_text: str


def run_propka_analysis(pdb_text: str, filename: str) -> PropkaAnalysis:
    molecule = propka.run.single(
        filename,
        stream=StringIO(pdb_text),
        write_pka=False,
    )
    conformation_name = (
        "AVR"
        if "AVR" in molecule.conformations
        else next(iter(molecule.conformations))
    )
    conformation = molecule.conformations[conformation_name]
    titration_sites = tuple(_extract_sites(conformation.get_titratable_groups()))
    charge_profile = tuple(
        ChargePoint(ph=ph, unfolded_charge=unfolded, folded_charge=folded)
        for ph, unfolded, folded in molecule.get_charge_profile(
            conformation=conformation_name,
            grid=(0.0, 14.0, 0.25),
        )
    )
    folding_profile_data, optimum, _, _ = molecule.get_folding_profile(
        conformation=conformation_name,
        reference="neutral",
        grid=(0.0, 14.0, 0.25),
    )
    folding_profile = tuple(
        FoldingEnergyPoint(ph=ph, folding_free_energy=energy)
        for ph, energy in (folding_profile_data or [])
    )
    optimum_ph: float | None = None
    optimum_energy: float | None = None
    if optimum is not None:
        optimum_ph, optimum_energy = optimum
    summary_text = propka.output.get_summary_section(
        molecule,
        conformation_name,
        conformation.parameters,
    )
    determinants_text = propka.output.get_determinant_section(
        molecule,
        conformation_name,
        conformation.parameters,
    )
    return PropkaAnalysis(
        source_name=filename,
        structure_name=Path(filename).stem,
        pdb_text=pdb_text,
        titration_sites=titration_sites,
        charge_profile=charge_profile,
        folding_profile=folding_profile,
        folding_optimum_ph=optimum_ph,
        folding_optimum_energy=optimum_energy,
        summary_text=summary_text,
        determinants_text=determinants_text,
    )


def _extract_sites(groups: Iterable[object]) -> list[TitrationSite]:
    sites: list[TitrationSite] = []
    for group in groups:
        insertion_code = getattr(group.atom, "icode", " ").strip()
        sites.append(
            TitrationSite(
                label=group.label.strip(),
                residue_type=group.residue_type.strip(),
                chain_id=group.atom.chain_id.strip() or "?",
                residue_number=int(group.atom.res_num),
                insertion_code=insertion_code,
                pka=float(group.pka_value),
                model_pka=float(group.model_pka),
                charged_state_charge=float(group.charge),
                buried_fraction=float(group.buried),
                interactions=_extract_interactions(group),
            )
        )
    return sorted(
        sites,
        key=lambda site: (
            site.chain_id,
            site.residue_number,
            site.insertion_code,
            site.label,
        ),
    )


def _extract_interactions(group: object) -> tuple[SiteInteraction, ...]:
    interactions: list[SiteInteraction] = []
    site_label = group.label.strip()
    determinants = getattr(group, "determinants", {})
    for interaction_kind, entries in determinants.items():
        for determinant in entries:
            partner = getattr(determinant, "group", None)
            base_partner = getattr(partner, "group", partner)
            if base_partner is None:
                continue
            partner_label = base_partner.label.strip()
            if partner_label == site_label:
                continue
            interactions.append(
                SiteInteraction(
                    partner_label=partner_label,
                    partner_residue_type=getattr(
                        base_partner, "residue_type", getattr(base_partner, "res_name", "")
                    ).strip(),
                    chain_id=base_partner.atom.chain_id.strip() or "?",
                    residue_number=int(base_partner.atom.res_num),
                    insertion_code=getattr(base_partner.atom, "icode", " ").strip(),
                    interaction_kind=interaction_kind,
                    contribution=float(determinant.value),
                )
            )

    return tuple(
        sorted(
            interactions,
            key=lambda interaction: (
                -abs(interaction.contribution),
                interaction.chain_id,
                interaction.residue_number,
                interaction.insertion_code,
                interaction.partner_label,
                interaction.interaction_kind,
            ),
        )
    )
