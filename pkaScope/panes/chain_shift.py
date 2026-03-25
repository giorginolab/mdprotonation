from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

from ..app_state import AppState
from ..chain_shift import (
    ChainPkaShift,
    compare_chain_pkas,
    compare_site_sets_pkas,
    create_chain_shift_plot_figure,
    extract_chains_pdb,
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
CHAIN_COMPARE_SCHEMA_VERSION = "2026-03-25-chain-compare-v2"
ADVANCED_LAST_RUN_KEY = "chain-shift-advanced-last-run"
ADVANCED_LAST_RESULT_KEY = "chain-shift-advanced-last-result"


@st.cache_data(show_spinner=False)
def cached_chain_set_analysis(
    pdb_text: str,
    source_name: str,
    chain_ids: tuple[str, ...],
    set_label: str,
    schema_version: str,
) -> PropkaAnalysis:
    del schema_version
    chain_pdb_text = extract_chains_pdb(pdb_text, chain_ids)
    chain_token = "-".join(_chain_id_token(chain_id) for chain_id in chain_ids)
    chain_source_name = f"{Path(source_name).stem}-{set_label}-{chain_token}.pdb"
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

    advanced_mode = st.checkbox(
        "Enable advanced chain selection",
        value=False,
        key="chain-shift-advanced-mode",
    )

    if advanced_mode:
        st.caption(
            "Select chain sets for the complex and apo states. "
            "Apo options are restricted to chains selected for complex."
        )
        complex_selection = st.pills(
            "Complex chains",
            options=list(chain_ids),
            default=list(chain_ids),
            selection_mode="multi",
            format_func=_chain_display_label,
            key="chain-shift-advanced-complex-pills",
        ) or []
        complex_chain_ids = tuple(
            chain_id
            for chain_id in chain_ids
            if chain_id in set(complex_selection)
        )
        apo_default = list(complex_chain_ids[:-1]) if len(complex_chain_ids) > 2 else []
        apo_selection = st.pills(
            "Apo chains",
            options=list(complex_chain_ids),
            default=apo_default,
            selection_mode="multi",
            format_func=_chain_display_label,
            key="chain-shift-advanced-apo-pills",
        ) or []
        apo_chain_ids = tuple(
            chain_id
            for chain_id in complex_chain_ids
            if chain_id in set(apo_selection)
        )
        st.caption(f"Complex set: `{_chain_set_label(complex_chain_ids)}`")
        st.caption(f"Apo set: `{_chain_set_label(apo_chain_ids)}`")
        if len(complex_chain_ids) <= 1:
            st.warning("Select at least 2 chains in Complex chains.")
            return
        if len(apo_chain_ids) <= 0:
            st.warning("Select at least 1 chain in Apo chains.")
            return
        if len(apo_chain_ids) >= len(complex_chain_ids):
            st.warning(
                "Apo chains must be a proper subset of Complex chains."
            )
            return

        run_clicked = st.button(
            "Run advanced comparison",
            type="primary",
            key="chain-shift-advanced-run-button",
        )
        if run_clicked:
            try:
                with st.spinner("Running PROPKA for selected chain sets..."):
                    complex_analysis = cached_chain_set_analysis(
                        analysis.pdb_text,
                        analysis.source_name,
                        complex_chain_ids,
                        "complex",
                        CHAIN_COMPARE_SCHEMA_VERSION,
                    )
                    monomer_analysis = cached_chain_set_analysis(
                        analysis.pdb_text,
                        analysis.source_name,
                        apo_chain_ids,
                        "apo",
                        CHAIN_COMPARE_SCHEMA_VERSION,
                    )
            except Exception as exc:
                st.error(
                    "Chain-set PROPKA failed for the advanced selection. "
                    "The comparison view could not be assembled."
                )
                st.exception(exc)
                return
            comparison = compare_site_sets_pkas(
                comparison_label="Complex vs Apo",
                complex_sites=complex_analysis.titration_sites,
                monomer_sites=monomer_analysis.titration_sites,
            )
            comparison_scope_key = (
                f"advanced-{_chain_set_key(complex_chain_ids)}-vs-{_chain_set_key(apo_chain_ids)}"
            )
            st.session_state[ADVANCED_LAST_RUN_KEY] = {
                "chain_ids": chain_ids,
                "complex_chain_ids": complex_chain_ids,
                "apo_chain_ids": apo_chain_ids,
            }
            st.session_state[ADVANCED_LAST_RESULT_KEY] = {
                "comparison": comparison,
                "comparison_scope_key": comparison_scope_key,
            }

        last_run = st.session_state.get(ADVANCED_LAST_RUN_KEY)
        last_result = st.session_state.get(ADVANCED_LAST_RESULT_KEY)
        if last_run is None or tuple(last_run.get("chain_ids", ())) != chain_ids:
            if ADVANCED_LAST_RUN_KEY in st.session_state:
                del st.session_state[ADVANCED_LAST_RUN_KEY]
            if ADVANCED_LAST_RESULT_KEY in st.session_state:
                del st.session_state[ADVANCED_LAST_RESULT_KEY]
            st.info(
                "Choose complex/apo chain sets and click `Run advanced comparison` "
                "to compute pKa deltas."
            )
            return

        if last_result is None:
            st.info(
                "Choose complex/apo chain sets and click `Run advanced comparison` "
                "to compute pKa deltas."
            )
            return

        run_complex_chain_ids = tuple(last_run["complex_chain_ids"])
        run_apo_chain_ids = tuple(last_run["apo_chain_ids"])
        if (
            run_complex_chain_ids != complex_chain_ids
            or run_apo_chain_ids != apo_chain_ids
        ):
            st.info(
                "Selections changed. Click Run advanced comparison to update results."
            )
            return

        comparison = last_result["comparison"]
        comparison_scope_key = last_result["comparison_scope_key"]
    else:
        selected_chain = st.selectbox(
            "Select chain",
            options=chain_ids,
            key="chain-shift-chain-select",
        )

        try:
            with st.spinner(f"Running chain-only PROPKA for chain {selected_chain}..."):
                monomer_analysis = cached_chain_set_analysis(
                    analysis.pdb_text,
                    analysis.source_name,
                    (selected_chain,),
                    "monomer",
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
        comparison_scope_key = selected_chain

    if not comparison.shifts:
        st.info(
            "No matched titratable sites were found between "
            "the selected complex and monomer calculations."
        )
        return

    matched_sites_col, max_shift_col, hot_spot_col = st.columns(3)
    matched_sites_col.metric("Matched titratable sites", f"{len(comparison.shifts)}")
    max_shift_col.metric("Max |Delta pKa|", f"{comparison.max_absolute_delta:.2f}")
    hot_spot_residue_count = len(top_shift_residue_keys(comparison.shifts, limit=12))
    hot_spot_col.metric("Hot-spot residues (|Delta| >= 0.5)", f"{hot_spot_residue_count}")

    st.plotly_chart(
        create_chain_shift_plot_figure(comparison.shifts, current_ph=app_state.ph),
        use_container_width=True,
        config={"displaylogo": False},
    )

    ordered_shifts = list(comparison.shifts)
    shift_df = _build_shift_dataframe(ordered_shifts)
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
    table_key = f"chain-shift-table-{comparison_scope_key}"
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
    if selected_row is None or selected_row >= len(ordered_shifts):
        selected_shift = ordered_shifts[0]
    else:
        selected_shift = ordered_shifts[selected_row]

    focus_targets = cached_chain_focus_targets(analysis.pdb_text)
    focus_residue = _resolve_focus_target(focus_targets, selected_shift)

    st.subheader("3D Shift Hot Spots")
    st.caption(
        "The viewer loads the full complex and automatically focuses on "
        "the first listed site. Select a table row to inspect another hot spot."
    )
    st_molstar_focusable_content(
        analysis.pdb_text,
        "pdb",
        file_name=f"{app_state.analysis.structure_name}-complex.pdb",
        focus_residue=focus_residue,
        height=f"{SHIFT_VIEWER_HEIGHT}px",
        key=f"chain-shift-viewer-{app_state.analysis.structure_name}-{comparison_scope_key}",
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
                "Site": _site_column_label(shift),
                "PROPKA label": shift.label,
                "Complex pKa": shift.complex_pka,
                "Monomer pKa": shift.monomer_pka,
                "Delta pKa": shift.delta_pka,
                "|Delta pKa|": shift.absolute_delta_pka,
            }
            for shift in shifts
        ]
    )


def _site_column_label(shift: ChainPkaShift) -> str:
    insertion_code = shift.insertion_code.strip()
    residue_token = (
        f"{shift.residue_number:>5}{insertion_code}"
        if insertion_code
        else f"{shift.residue_number:>5}"
    )
    return f"{shift.chain_id}:{residue_token}-{shift.residue_type}"


def _chain_set_label(chain_ids: tuple[str, ...]) -> str:
    if not chain_ids:
        return "(none)"
    return ", ".join(_chain_display_label(chain_id) for chain_id in chain_ids)


def _chain_set_key(chain_ids: tuple[str, ...]) -> str:
    if not chain_ids:
        return "none"
    return "-".join(_chain_id_token(chain_id) for chain_id in chain_ids)


def _chain_display_label(chain_id: str) -> str:
    return chain_id if chain_id != "?" else "(blank)"


def _chain_id_token(chain_id: str) -> str:
    return "blank" if chain_id == "?" else chain_id


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
