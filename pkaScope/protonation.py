from __future__ import annotations

from dataclasses import dataclass
from math import exp
from statistics import mean

from .pdb_utils import atom_residue_key, load_pdb_universe, write_pdb_atoms
from .propka_analysis import TitrationSite


@dataclass(frozen=True)
class SiteState:
    site: TitrationSite
    ph: float
    protonated_fraction: float
    deprotonated_fraction: float
    current_charge: float
    transition_score: float
    dominant_state: str


@dataclass(frozen=True)
class ResidueVisualEncoding:
    chain_id: str
    residue_number: int
    insertion_code: str
    average_protonated_fraction: float
    transition_intensity: float
    net_charge: float
    site_count: int

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain_id, self.residue_number, self.insertion_code)


def protonated_fraction(ph: float, pka: float) -> float:
    exponent = (ph - pka) * 2.302585092994046
    if exponent > 700:
        return 0.0
    if exponent < -700:
        return 1.0
    return 1.0 / (1.0 + exp(exponent))


def evaluate_sites(sites: tuple[TitrationSite, ...], ph: float) -> list[SiteState]:
    states: list[SiteState] = []
    for site in sites:
        protonated = protonated_fraction(ph, site.pka)
        deprotonated = 1.0 - protonated
        current_charge = _charge_at_fraction(site, protonated)
        transition_score = 1.0 - abs((protonated * 2.0) - 1.0)
        states.append(
            SiteState(
                site=site,
                ph=ph,
                protonated_fraction=protonated,
                deprotonated_fraction=deprotonated,
                current_charge=current_charge,
                transition_score=transition_score,
                dominant_state=_dominant_state(site, protonated),
            )
        )
    return states


def site_charge_at_ph(site: TitrationSite, ph: float) -> float:
    protonated = protonated_fraction(ph, site.pka)
    return _charge_at_fraction(site, protonated)


def summarize_residues(site_states: list[SiteState]) -> dict[tuple[str, int, str], ResidueVisualEncoding]:
    grouped: dict[tuple[str, int, str], list[SiteState]] = {}
    for state in site_states:
        grouped.setdefault(state.site.residue_key, []).append(state)

    return {
        residue_key: ResidueVisualEncoding(
            chain_id=residue_key[0],
            residue_number=residue_key[1],
            insertion_code=residue_key[2],
            average_protonated_fraction=mean(
                state.protonated_fraction for state in residue_states
            ),
            transition_intensity=max(
                state.transition_score for state in residue_states
            ),
            net_charge=sum(state.current_charge for state in residue_states),
            site_count=len(residue_states),
        )
        for residue_key, residue_states in grouped.items()
    }


def render_ph_encoded_pdb(
    pdb_text: str,
    residue_encodings: dict[tuple[str, int, str], ResidueVisualEncoding],
) -> str:
    universe = load_pdb_universe(pdb_text)
    for atom in universe.atoms:
        residue_key = atom_residue_key(atom)
        encoding = residue_encodings.get(residue_key)
        if encoding is None:
            continue
        occupancy = min(max(encoding.average_protonated_fraction, 0.0), 1.0)
        b_factor = min(max(encoding.transition_intensity * 100.0, 0.0), 100.0)
        atom.occupancy = occupancy
        atom.tempfactor = b_factor

    encoded_lines = [
        "REMARK 950 OCCUPANCY STORES AVERAGE PROTONATED FRACTION (0.00-1.00)\n",
        "REMARK 950 B-FACTOR STORES TRANSITION INTENSITY * 100 (0.00-100.00)\n",
        write_pdb_atoms(universe.atoms),
    ]
    return "".join(encoded_lines)


def _charge_at_fraction(site: TitrationSite, protonated: float) -> float:
    if site.charged_state_charge > 0:
        return protonated * site.charged_state_charge
    return (1.0 - protonated) * site.charged_state_charge


def _dominant_state(site: TitrationSite, protonated: float) -> str:
    if 0.1 < protonated < 0.9:
        return "Transitioning"
    if site.charged_state_charge > 0:
        return "Cationic" if protonated >= 0.5 else "Neutral"
    return "Neutral" if protonated >= 0.5 else "Anionic"
