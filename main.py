from __future__ import annotations

from functools import partial
from pathlib import Path

import pandas as pd
import streamlit as st

from mdprotonation.propka_analysis import PropkaAnalysis, run_propka_analysis
from mdprotonation.protonation import (
    evaluate_sites,
    render_ph_encoded_pdb,
    summarize_residues,
)
from mdprotonation.viewer import (
    compute_residue_focus_targets,
    st_molstar_focusable_content,
)


st.set_page_config(
    page_title="Protein Protonation Explorer",
    layout="wide",
)

ANALYSIS_SCHEMA_VERSION = "2026-03-23-profiles-v1"


@st.cache_data(show_spinner=False)
def cached_propka_analysis(
    pdb_text: str,
    filename: str,
    schema_version: str,
) -> PropkaAnalysis:
    del schema_version
    return run_propka_analysis(pdb_text, filename)


def load_propka_analysis(pdb_text: str, filename: str) -> PropkaAnalysis:
    analysis = cached_propka_analysis(
        pdb_text,
        filename,
        ANALYSIS_SCHEMA_VERSION,
    )
    required_fields = (
        "folding_profile",
        "folding_optimum_ph",
        "folding_optimum_energy",
    )
    if all(hasattr(analysis, field) for field in required_fields):
        return analysis

    cached_propka_analysis.clear()
    return cached_propka_analysis(
        pdb_text,
        filename,
        ANALYSIS_SCHEMA_VERSION,
    )


def style_site_table(row: pd.Series) -> list[str]:
    charge = float(row["Charge"])

    # Inverted five-bin scale across charge:
    # red -> pink -> white -> slate blue -> deep blue
    if charge <= -0.75:
        background = "#c62828"  # red
        foreground = "#ffffff"
    elif charge <= -0.25:
        background = "#f6c1d6"  # pink
        foreground = "#3f1025"
    elif charge < 0.25:
        background = "#ffffff"  # white
        foreground = "#111827"
    elif charge < 0.75:
        background = "#6f82be"  # slate blue
        foreground = "#ffffff"
    else:
        background = "#0d47a1"  # deep blue
        foreground = "#ffffff"

    return [f"background-color: {background}; color: {foreground}"] * len(row)


def selected_rows_from_selection(selection: object | None) -> list[int]:
    if selection is None:
        return []

    rows = getattr(selection, "rows", None)
    if rows is None and isinstance(selection, dict):
        rows = selection.get("rows")
    if isinstance(rows, list):
        parsed_rows = [row for row in rows if isinstance(row, int)]
        if parsed_rows:
            return parsed_rows

    cells = getattr(selection, "cells", None)
    if cells is None and isinstance(selection, dict):
        cells = selection.get("cells")
    if not isinstance(cells, list):
        return []

    cell_rows: list[int] = []
    for cell in cells:
        if isinstance(cell, (tuple, list)) and cell and isinstance(cell[0], int):
            cell_rows.append(cell[0])

    unique_rows: list[int] = []
    for row in cell_rows:
        if row not in unique_rows:
            unique_rows.append(row)
    return unique_rows


def selected_rows_from_state(key: str) -> list[int]:
    state = st.session_state.get(key)
    if state is None:
        return []

    selection = getattr(state, "selection", None)
    if selection is None and isinstance(state, dict):
        selection = state.get("selection")
    return selected_rows_from_selection(selection)


def set_active_selection_source(source: str) -> None:
    st.session_state["active_selection_source"] = source


def main() -> None:
    st.title("Protein Protonation Explorer")
    st.caption(
        "Explore continuous residue protonation as a function of pH using "
        "[PROPka](https://github.com/jensengroup/propka)."
    )

    with st.sidebar:
        st.header("Structure")
        source_mode = st.radio(
            "Input source",
            options=("Example PDB", "Upload PDB"),
            index=0,
        )

        if source_mode == "Example PDB":
            example_paths = sorted(Path("examples").glob("*.pdb"))
            if not example_paths:
                st.error("No example PDB files were found in examples/.")
                return
            selected_example = st.selectbox(
                "Example structure",
                options=example_paths,
                format_func=lambda path: path.name,
            )
            source_name = selected_example.name
            pdb_text = selected_example.read_text(encoding="utf-8")
        else:
            upload = st.file_uploader("Upload a PDB file", type=["pdb"])
            if upload is None:
                st.info("Upload a PDB file to run PROPka.")
                return
            source_name = upload.name
            pdb_text = upload.getvalue().decode("utf-8")

        st.header("pH")
        ph = st.slider(
            "Solution pH",
            min_value=0.0,
            max_value=14.0,
            value=7.0,
            step=0.1,
        )
        transition_only = st.checkbox("Show only transitioning sites", value=False)

    try:
        with st.spinner("Running PROPka and assembling titration states..."):
            analysis = load_propka_analysis(pdb_text, source_name)
    except Exception as exc:
        st.error("PROPka could not process this structure.")
        st.exception(exc)
        return

    site_states = evaluate_sites(analysis.titration_sites, ph)
    standard_site_states = evaluate_sites(analysis.titration_sites, 7.0)
    standard_state_by_site = {
        state.site.site_key: state.dominant_state
        for state in standard_site_states
    }
    residue_encodings = summarize_residues(site_states)
    encoded_pdb_text = render_ph_encoded_pdb(analysis.pdb_text, residue_encodings)
    residue_focus_targets = compute_residue_focus_targets(analysis.pdb_text)

    transitioning_count = sum(
        1 for state in site_states if state.dominant_state == "Transitioning"
    )
    folded_charge = sum(state.current_charge for state in site_states)
    responsive_sites = sorted(
        site_states,
        key=lambda state: (abs(state.ph - state.site.pka), -state.transition_score),
    )
    top_responsive_sites = responsive_sites[:8]
    responsive_selected_state = None

    with st.sidebar:
        st.subheader("Current State")
        metric_col1, metric_col2 = st.columns(2)
        metric_col1.metric("Titratable sites", len(site_states))
        metric_col2.metric("Transitioning", transitioning_count)
        metric_col1.metric("Folded charge", f"{folded_charge:.2f}")
        metric_col2.metric("Structure", analysis.structure_name)

        st.markdown("#### Residues responding most at this pH")
        responsive_rows = [
            {
                "Site": state.site.label,
                "pKa": state.site.pka,
                "State": state.dominant_state,
                "% Protonated": state.protonated_fraction * 100.0,
            }
            for state in top_responsive_sites
        ]
        responsive_df = pd.DataFrame(responsive_rows)
        responsive_table_event = st.dataframe(
            responsive_df,
            use_container_width=True,
            hide_index=True,
            key="responsive-site-table",
            on_select=partial(set_active_selection_source, "responsive"),
            selection_mode="single-cell",
            column_config={
                "pKa": st.column_config.NumberColumn(format="%.2f"),
                "% Protonated": st.column_config.NumberColumn(format="%.0f%%"),
            },
        )
        responsive_selected_rows = (
            selected_rows_from_selection(responsive_table_event.selection)
            or selected_rows_from_state("responsive-site-table")
        )
        if (
            responsive_selected_rows
            and responsive_selected_rows[0] < len(top_responsive_sites)
        ):
            responsive_selected_state = top_responsive_sites[responsive_selected_rows[0]]

    overview_tab, profiles_tab, propka_tab = st.tabs(
        ["Explorer", "Profiles", "PROPka Data"]
    )

    with overview_tab:
        filtered_states = [
            state
            for state in site_states
            if not transition_only or state.dominant_state == "Transitioning"
        ]
        sorted_states = sorted(
            filtered_states,
            key=lambda state: (
                state.site.chain_id,
                state.site.residue_number,
                state.site.insertion_code,
                state.site.label,
            ),
        )
        selected_rows = selected_rows_from_state("site-state-table")
        main_table_selected_state = (
            sorted_states[selected_rows[0]]
            if selected_rows and selected_rows[0] < len(sorted_states)
            else None
        )
        active_selection_source = st.session_state.get("active_selection_source", "main")
        if (
            active_selection_source == "responsive"
            and responsive_selected_state is not None
        ):
            selected_state = responsive_selected_state
        elif main_table_selected_state is not None:
            selected_state = main_table_selected_state
        else:
            selected_state = responsive_selected_state

        focus_residue = (
            residue_focus_targets.get(selected_state.site.residue_key)
            if selected_state is not None
            else None
        )

        st.subheader("Structure Viewer")
        st_molstar_focusable_content(
            encoded_pdb_text,
            "pdb",
            file_name=f"{analysis.structure_name}-ph-{ph:.1f}.pdb",
            focus_residue=focus_residue,
            height="640px",
            key=f"molstar-{analysis.structure_name}-{ph:.1f}",
        )
        st.caption(
            "Viewer encoding: occupancy stores average protonated fraction per "
            "titratable residue and B-factor stores transition intensity x100."
        )
        if selected_state is not None:
            st.caption(
                f"Camera focus target: `{selected_state.site.label}` "
                "from the selected table row."
            )

        st.subheader("Residue Protonation States")
        state_rows = [
            {
                "pH7": (
                    ""
                    if state.dominant_state
                    == standard_state_by_site.get(state.site.site_key, state.dominant_state)
                    else "⚠️"
                ),
                "Site": state.site.label,
                "Residue": state.site.residue_type,
                "Chain": state.site.chain_id,
                "Residue #": state.site.residue_number,
                "pKa": state.site.pka,
                "Model pKa": state.site.model_pka,
                "Buried fraction": state.site.buried_fraction,
                "% Protonated": state.protonated_fraction * 100.0,
                "Charge": state.current_charge,
                "State": state.dominant_state,
            }
            for state in sorted_states
        ]
        state_df = pd.DataFrame(state_rows)
        styled_state_df = state_df.style.apply(style_site_table, axis=1)
        styled_state_df = styled_state_df.format(
            {
                "pKa": "{:.2f}",
                "Model pKa": "{:.2f}",
                "Buried fraction": "{:.3f}",
                "Charge": "{:.3f}",
                "% Protonated": "{:.0f}",
            }
        )
        table_event = st.dataframe(
            styled_state_df,
            use_container_width=True,
            hide_index=True,
            key="site-state-table",
            on_select=partial(set_active_selection_source, "main"),
            selection_mode="single-cell",
        )
        selected_rows = selected_rows_from_selection(table_event.selection) or selected_rows
        if selected_rows and selected_rows[0] < len(sorted_states):
            main_table_selected_state = sorted_states[selected_rows[0]]
            if st.session_state.get("active_selection_source", "main") == "main":
                selected_state = main_table_selected_state
        elif selected_state is None:
            selected_state = responsive_selected_state

        st.caption(
            "Row colors by charge (inverted scale): <= -0.75 red, "
            "(-0.75, -0.25] pink, (-0.25, +0.25) white, "
            "[+0.25, +0.75) slate blue, >= +0.75 deep blue."
        )
        st.caption("`⚠️` marks sites not in their pH 7 dominant state.")

        st.subheader("Selected Site Interactions")
        if selected_state is not None:
            st.write(
                (
                    f"`{selected_state.site.label}` | "
                    f"buried fraction `{selected_state.site.buried_fraction:.3f}` | "
                    f"pKa `{selected_state.site.pka:.2f}`"
                )
            )
            interaction_rows = [
                {
                    "Partner site": interaction.partner_label,
                    "Chain": interaction.chain_id,
                    "Residue #": interaction.residue_number,
                    "Kind": interaction.interaction_kind,
                    "Contribution": round(interaction.contribution, 3),
                }
                for interaction in selected_state.site.interactions
            ]
            if interaction_rows:
                st.dataframe(
                    interaction_rows,
                    use_container_width=True,
                    hide_index=True,
                )
            else:
                st.info("PROPka reported no non-self interaction residues for this site.")
        else:
            st.info(
                "Click a row in either residue table to inspect PROPka interaction residues."
            )

    with profiles_tab:
        st.subheader("Charge Profile")
        charge_points = [
            {
                "pH": point.ph,
                "Folded charge": point.folded_charge,
                "Unfolded charge": point.unfolded_charge,
            }
            for point in analysis.charge_profile
        ]
        st.line_chart(charge_points, x="pH", y=["Folded charge", "Unfolded charge"])

        st.subheader("Free Energy of Folding")
        if analysis.folding_profile:
            if (
                analysis.folding_optimum_ph is not None
                and analysis.folding_optimum_energy is not None
            ):
                opt_col1, opt_col2 = st.columns(2)
                opt_col1.metric("Optimum pH", f"{analysis.folding_optimum_ph:.2f}")
                opt_col2.metric(
                    "Minimum folding free energy",
                    f"{analysis.folding_optimum_energy:.2f} kcal/mol",
                )

            folding_points = [
                {
                    "pH": point.ph,
                    "Free energy of folding (kcal/mol)": point.folding_free_energy,
                }
                for point in analysis.folding_profile
            ]
            st.line_chart(
                folding_points,
                x="pH",
                y=["Free energy of folding (kcal/mol)"],
            )
        else:
            st.info("PROPka did not return a folding free-energy profile for this structure.")

    with propka_tab:
        st.subheader("PROPka Summary")
        st.code(analysis.summary_text.strip(), language="text")
        st.subheader("PROPka Determinants")
        st.code(analysis.determinants_text.strip(), language="text")


if __name__ == "__main__":
    main()
