from __future__ import annotations

from pathlib import Path

import streamlit as st

from mdprotonation.app_state import build_app_state
from mdprotonation.online_pdb import OnlinePdbError, fetch_pdb_from_rcsb, normalize_pdb_id
from mdprotonation.panes import (
    render_chain_shift_tab,
    render_explorer_tab,
    render_pka_plot_tab,
    render_profiles_tab,
    render_propka_data_tab,
    render_sidebar_summary,
)
from mdprotonation.propka_analysis import PropkaAnalysis, run_propka_analysis


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


@st.cache_data(show_spinner=False)
def cached_online_pdb_fetch(pdb_id: str) -> str:
    return fetch_pdb_from_rcsb(pdb_id)


def main() -> None:
    st.title("Protein Protonation Explorer")
    st.caption(
        "Explore continuous residue protonation as a function of pH using "
        "[PROPKA](https://github.com/jensengroup/propka). "
        "[giorginolab.it](https://www.giorginolab.it)"
    )

    with st.sidebar:
        st.header("Structure")
        source_mode = st.radio(
            "Input source",
            options=("Example PDB", "Upload PDB", "PDB Online"),
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
        elif source_mode == "Upload PDB":
            upload = st.file_uploader("Upload a PDB file", type=["pdb"])
            if upload is None:
                st.info("Upload a PDB file to run PROPKA.")
                return
            source_name = upload.name
            pdb_text = upload.getvalue().decode("utf-8")
        else:
            pdb_id_input = st.text_input("RCSB PDB ID", value="7BCQ", max_chars=4)
            if not pdb_id_input.strip():
                st.info("Enter a 4-character PDB ID to download a structure.")
                return
            try:
                pdb_id = normalize_pdb_id(pdb_id_input)
            except ValueError as exc:
                st.warning(str(exc))
                return
            try:
                with st.spinner(f"Downloading {pdb_id} from RCSB..."):
                    pdb_text = cached_online_pdb_fetch(pdb_id)
            except OnlinePdbError as exc:
                st.error(str(exc))
                return
            source_name = f"{pdb_id}.pdb"

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
        with st.spinner("Running PROPKA and assembling titration states..."):
            analysis = load_propka_analysis(pdb_text, source_name)
    except Exception as exc:
        st.error("PROPKA could not process this structure.")
        st.exception(exc)
        return

    app_state = build_app_state(analysis, ph)
    responsive_selected_state = render_sidebar_summary(app_state)

    overview_tab, pka_plot_tab, profiles_tab, chain_shift_tab, propka_tab = st.tabs(
        ["Explorer", "pKa Plot", "Profiles", "Chain Shifts", "PROPKA Data"]
    )

    with overview_tab:
        render_explorer_tab(
            app_state,
            transition_only=transition_only,
            responsive_selected_state=responsive_selected_state,
        )

    with pka_plot_tab:
        render_pka_plot_tab(app_state)

    with profiles_tab:
        render_profiles_tab(analysis)

    with chain_shift_tab:
        render_chain_shift_tab(analysis, app_state)

    with propka_tab:
        render_propka_data_tab(analysis)

    st.divider()
    st.markdown(
        "Created by [Toni Giorgino](https://www.giorginolab.it) | AI coded | "
        "[GitHub](https://github.com/giorginolab/pkaScope#) | "
        "[Hugging Face Space](https://huggingface.co/spaces/tonigi/pkaScope)"
    )


if __name__ == "__main__":
    main()
