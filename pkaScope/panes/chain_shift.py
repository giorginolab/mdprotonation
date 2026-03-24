from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from ..app_state import AppState
from ..chain_shift import (
    ChainPkaShift,
    compare_chain_pkas,
    create_chain_shift_plot_figure,
    extract_chain_pdb,
    list_chain_ids,
    top_shift_residue_keys,
)
from ..propka_analysis import PropkaAnalysis, run_propka_analysis
from ..ui_state import get_table_selection
from ..viewer import (
    ResidueFocusTarget,
    compute_residue_focus_targets,
    st_molstar_focusable_content,
)

SHIFT_VIEWER_HEIGHT = 620
SHIFT_TABLE_HEIGHT = 380
CHAIN_COMPARE_SCHEMA_VERSION = "2026-03-24-chain-compare-v1"


@st.cache_data(show_spinner=False)
def cached_chain_monomer_analysis(
    pdb_text: str,
    source_name: str,
    chain_id: str,
    schema_version: str,
) -> PropkaAnalysis:
    del schema_version
    chain_pdb_text = extract_chain_pdb(pdb_text, chain_id)
    chain_source_name = f"{Path(source_name).stem}-chain-{chain_id}.pdb"
    return run_propka_analysis(chain_pdb_text, chain_source_name)


@st.cache_data(show_spinner=False)
def cached_chain_focus_targets(
    pdb_text: str,
) -> dict[tuple[str, int, str], ResidueFocusTarget]:
    return compute_residue_focus_targets(pdb_text)


def render_chain_shift_tab(
    analysis: PropkaAnalysis,
    app_state: AppState,
) -> None:
    st.subheader("Chain pKa Shifts (Complex vs Monomer)")
    chain_ids = list_chain_ids(analysis.pdb_text)
    if not chain_ids:
        st.info("No chain IDs were found in this structure.")
        return

    selected_chain = st.selectbox(
        "Select chain",
        options=chain_ids,
        key="chain-shift-chain-select",
    )

    try:
        with st.spinner(f"Running chain-only PROPKA for chain {selected_chain}..."):
            monomer_analysis = cached_chain_monomer_analysis(
                analysis.pdb_text,
                analysis.source_name,
                selected_chain,
                CHAIN_COMPARE_SCHEMA_VERSION,
            )
    except Exception as exc:
        st.error(
            f"Chain-only PROPKA failed for chain {selected_chain}. "
            "The comparison view could not be assembled."
        )
        st.exception(exc)
        return

    comparison = compare_chain_pkas(
        chain_id=selected_chain,
        complex_sites=analysis.titration_sites,
        monomer_sites=monomer_analysis.titration_sites,
    )
    if not comparison.shifts:
        st.info(
            "No matched titratable sites were found for the selected chain between "
            "complex and monomer calculations."
        )
        return

    matched_sites_col, max_shift_col, hot_spot_col = st.columns(3)
    matched_sites_col.metric("Matched titratable sites", f"{len(comparison.shifts)}")
    max_shift_col.metric("Max |Delta pKa|", f"{comparison.max_absolute_delta:.2f}")
    hot_spot_residue_count = len(top_shift_residue_keys(comparison.shifts, limit=12))
    hot_spot_col.metric("Hot-spot residues (|Delta| >= 0.5)", f"{hot_spot_residue_count}")

    st.pyplot(
        create_chain_shift_plot_figure(comparison.shifts, current_ph=app_state.ph),
        use_container_width=True,
    )

    ranked_shifts = sorted(
        comparison.shifts,
        key=lambda shift: (
            -shift.absolute_delta_pka,
            shift.residue_number,
            shift.insertion_code,
            shift.label,
        ),
    )
    shift_df = _build_shift_dataframe(ranked_shifts)
    styled_shift_df = shift_df.style.format(
        {
            "Complex pKa": "{:.2f}",
            "Monomer pKa": "{:.2f}",
            "Delta pKa": "{:+.2f}",
            "|Delta pKa|": "{:.2f}",
        }
    ).background_gradient(
        subset=["|Delta pKa|"],
        cmap="Reds",
    )
    table_key = f"chain-shift-table-{selected_chain}"
    table_event = st.dataframe(
        styled_shift_df,
        hide_index=True,
        use_container_width=True,
        height=SHIFT_TABLE_HEIGHT,
        key=table_key,
        selection_mode="single-row",
        on_select="rerun",
    )

    selected_row = get_table_selection(table_key, table_event.selection).first_row
    if selected_row is None or selected_row >= len(ranked_shifts):
        selected_shift = ranked_shifts[0]
    else:
        selected_shift = ranked_shifts[selected_row]

    focus_targets = cached_chain_focus_targets(analysis.pdb_text)
    focus_residue = _resolve_focus_target(focus_targets, selected_shift)

    st.subheader("3D Shift Hot Spots")
    st.caption(
        "The viewer loads the full complex and automatically focuses on "
        "the highest-shift site. Select a table row to inspect another hot spot."
    )
    st_molstar_focusable_content(
        analysis.pdb_text,
        "pdb",
        file_name=f"{app_state.analysis.structure_name}-complex.pdb",
        focus_residue=focus_residue,
        height=f"{SHIFT_VIEWER_HEIGHT}px",
        key=f"chain-shift-viewer-{app_state.analysis.structure_name}-{selected_chain}",
    )
    st.caption(
        f"Focused site: `{selected_shift.residue_label}` "
        f"(Delta pKa {selected_shift.delta_pka:+.2f})."
    )

    if comparison.unmatched_complex_labels:
        st.caption(
            "Complex-only titration sites not matched in monomer: "
            + ", ".join(comparison.unmatched_complex_labels[:8])
            + (
                "..."
                if len(comparison.unmatched_complex_labels) > 8
                else ""
            )
        )
    if comparison.unmatched_monomer_labels:
        st.caption(
            "Monomer-only titration sites not matched in complex: "
            + ", ".join(comparison.unmatched_monomer_labels[:8])
            + (
                "..."
                if len(comparison.unmatched_monomer_labels) > 8
                else ""
            )
        )


def _build_shift_dataframe(shifts: list[ChainPkaShift]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Site": shift.residue_label,
                "PROPKA label": shift.label,
                "Complex pKa": shift.complex_pka,
                "Monomer pKa": shift.monomer_pka,
                "Delta pKa": shift.delta_pka,
                "|Delta pKa|": shift.absolute_delta_pka,
            }
            for shift in shifts
        ]
    )


def _resolve_focus_target(
    focus_targets: dict[tuple[str, int, str], ResidueFocusTarget],
    shift: ChainPkaShift,
) -> ResidueFocusTarget | None:
    focused = focus_targets.get(shift.residue_key)
    if focused is not None:
        return focused
    for residue_key, target in focus_targets.items():
        if residue_key[1] == shift.residue_number and residue_key[2] == shift.insertion_code:
            return target
    return None
