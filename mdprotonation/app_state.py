from __future__ import annotations

from dataclasses import dataclass

from .propka_analysis import PropkaAnalysis
from .protonation import (
    ResidueVisualEncoding,
    SiteState,
    evaluate_sites,
    render_ph_encoded_pdb,
    summarize_residues,
)
from .viewer import ResidueFocusTarget, compute_residue_focus_targets

SiteKey = tuple[str, int, str, str]


def site_state_sort_key(state: SiteState) -> tuple[str, int, str, str]:
    return (
        state.site.chain_id,
        state.site.residue_number,
        state.site.insertion_code,
        state.site.label,
    )


def responsive_site_sort_key(state: SiteState) -> tuple[float, float]:
    return (abs(state.ph - state.site.pka), -state.transition_score)


@dataclass(frozen=True)
class AppState:
    analysis: PropkaAnalysis
    ph: float
    site_states: tuple[SiteState, ...]
    sorted_site_states: tuple[SiteState, ...]
    standard_state_by_site: dict[SiteKey, str]
    residue_encodings: dict[tuple[str, int, str], ResidueVisualEncoding]
    encoded_pdb_text: str
    residue_focus_targets: dict[tuple[str, int, str], ResidueFocusTarget]
    transitioning_count: int
    folded_charge: float
    top_responsive_sites: tuple[SiteState, ...]

    def filtered_site_states(self, transition_only: bool) -> tuple[SiteState, ...]:
        if not transition_only:
            return self.sorted_site_states
        return tuple(
            state
            for state in self.sorted_site_states
            if state.dominant_state == "Transitioning"
        )


def build_app_state(analysis: PropkaAnalysis, ph: float) -> AppState:
    site_states = tuple(evaluate_sites(analysis.titration_sites, ph))
    sorted_site_states = tuple(sorted(site_states, key=site_state_sort_key))
    standard_state_by_site = {
        state.site.site_key: state.dominant_state
        for state in evaluate_sites(analysis.titration_sites, 7.0)
    }
    residue_encodings = summarize_residues(list(site_states))
    encoded_pdb_text = render_ph_encoded_pdb(analysis.pdb_text, residue_encodings)
    residue_focus_targets = compute_residue_focus_targets(analysis.pdb_text)
    transitioning_count = sum(
        1 for state in site_states if state.dominant_state == "Transitioning"
    )
    folded_charge = sum(state.current_charge for state in site_states)
    top_responsive_sites = tuple(
        sorted(site_states, key=responsive_site_sort_key)[:8]
    )

    return AppState(
        analysis=analysis,
        ph=ph,
        site_states=site_states,
        sorted_site_states=sorted_site_states,
        standard_state_by_site=standard_state_by_site,
        residue_encodings=residue_encodings,
        encoded_pdb_text=encoded_pdb_text,
        residue_focus_targets=residue_focus_targets,
        transitioning_count=transitioning_count,
        folded_charge=folded_charge,
        top_responsive_sites=top_responsive_sites,
    )
