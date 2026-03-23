from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Iterable

import propka.output
import propka.run


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
class PropkaAnalysis:
    source_name: str
    structure_name: str
    pdb_text: str
    titration_sites: tuple[TitrationSite, ...]
    charge_profile: tuple[ChargePoint, ...]
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
