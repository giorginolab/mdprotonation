from __future__ import annotations

from pathlib import Path

import streamlit as st
from streamlit_molstar import st_molstar_content

from mdprotonation.propka_analysis import PropkaAnalysis, run_propka_analysis
from mdprotonation.protonation import (
    evaluate_sites,
    render_ph_encoded_pdb,
    summarize_residues,
)


st.set_page_config(
    page_title="Protein Protonation Explorer",
    layout="wide",
)


@st.cache_data(show_spinner=False)
def cached_propka_analysis(pdb_text: str, filename: str) -> PropkaAnalysis:
    return run_propka_analysis(pdb_text, filename)


def main() -> None:
    st.title("Protein Protonation Explorer")
    st.caption(
        "Explore continuous residue protonation as a function of pH using PROPka."
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
            analysis = cached_propka_analysis(pdb_text, source_name)
    except Exception as exc:
        st.error("PROPka could not process this structure.")
        st.exception(exc)
        return

    site_states = evaluate_sites(analysis.titration_sites, ph)
    residue_encodings = summarize_residues(site_states)
    encoded_pdb_text = render_ph_encoded_pdb(analysis.pdb_text, residue_encodings)

    transitioning_count = sum(
        1 for state in site_states if state.dominant_state == "Transitioning"
    )
    folded_charge = sum(state.current_charge for state in site_states)
    responsive_sites = sorted(
        site_states,
        key=lambda state: (abs(state.ph - state.site.pka), -state.transition_score),
    )

    overview_tab, profiles_tab, propka_tab = st.tabs(
        ["Explorer", "Profiles", "PROPka Data"]
    )

    with overview_tab:
        viewer_col, detail_col = st.columns([1.75, 1.0], gap="large")

        with viewer_col:
            st.subheader("Structure Viewer")
            st_molstar_content(
                encoded_pdb_text,
                "pdb",
                file_name=f"{analysis.structure_name}-ph-{ph:.1f}.pdb",
                height="640px",
                key=f"molstar-{analysis.structure_name}-{ph:.1f}",
            )
            st.caption(
                "Viewer encoding: occupancy stores average protonated fraction per "
                "titratable residue and B-factor stores transition intensity x100."
            )

        with detail_col:
            st.subheader("Current State")
            metric_col1, metric_col2 = st.columns(2)
            metric_col1.metric("Titratable sites", len(site_states))
            metric_col2.metric("Transitioning", transitioning_count)
            metric_col1.metric("Folded charge", f"{folded_charge:.2f}")
            metric_col2.metric("Structure", analysis.structure_name)

            st.markdown("#### Residues responding most at this pH")
            for state in responsive_sites[:8]:
                st.write(
                    (
                        f"`{state.site.label}` | pKa {state.site.pka:.2f} | "
                        f"{state.dominant_state} | "
                        f"{state.protonated_fraction * 100:.0f}% protonated"
                    )
                )

        st.subheader("Residue Protonation States")
        filtered_states = [
            state
            for state in site_states
            if not transition_only or state.dominant_state == "Transitioning"
        ]
        state_rows = [
            {
                "Site": state.site.label,
                "Residue": state.site.residue_type,
                "Chain": state.site.chain_id,
                "Residue #": state.site.residue_number,
                "pKa": round(state.site.pka, 2),
                "Model pKa": round(state.site.model_pka, 2),
                "Buried fraction": round(state.site.buried_fraction, 3),
                "% Protonated": round(state.protonated_fraction * 100.0, 1),
                "Charge": round(state.current_charge, 3),
                "State": state.dominant_state,
            }
            for state in sorted(
                filtered_states,
                key=lambda state: (
                    state.site.chain_id,
                    state.site.residue_number,
                    state.site.insertion_code,
                    state.site.label,
                ),
            )
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
        table_event = st.dataframe(
            state_rows,
            use_container_width=True,
            hide_index=True,
            key="site-state-table",
            on_select="rerun",
            selection_mode="single-row",
        )

        selected_rows = table_event.selection.rows
        st.subheader("Selected Site Interactions")
        if selected_rows:
            selected_state = sorted_states[selected_rows[0]]
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
            st.info("Click a row in the site table to inspect PROPka interaction residues.")

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
