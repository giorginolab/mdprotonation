from __future__ import annotations

import streamlit as st

from ..propka_analysis import PropkaAnalysis


def render_profiles_tab(analysis: PropkaAnalysis) -> None:
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
    if not analysis.folding_profile:
        st.info("PROPKA did not return a folding free-energy profile for this structure.")
        return

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
