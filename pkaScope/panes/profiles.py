from __future__ import annotations

import plotly.graph_objects as go
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
    charge_figure = go.Figure()
    charge_figure.add_trace(
        go.Scatter(
            x=[point["pH"] for point in charge_points],
            y=[point["Folded charge"] for point in charge_points],
            mode="lines",
            name="Folded charge",
            line={"color": "#1d4ed8", "width": 2.2},
        )
    )
    charge_figure.add_trace(
        go.Scatter(
            x=[point["pH"] for point in charge_points],
            y=[point["Unfolded charge"] for point in charge_points],
            mode="lines",
            name="Unfolded charge",
            line={"color": "#ea580c", "width": 2.2},
        )
    )
    charge_figure.update_layout(
        margin={"l": 24, "r": 24, "t": 24, "b": 24},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        hovermode="x unified",
    )
    charge_figure.update_xaxes(title_text="pH")
    charge_figure.update_yaxes(title_text="Charge")
    st.plotly_chart(
        charge_figure,
        width="stretch",
        config={"displaylogo": False},
    )

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
    folding_figure = go.Figure()
    folding_figure.add_trace(
        go.Scatter(
            x=[point["pH"] for point in folding_points],
            y=[point["Free energy of folding (kcal/mol)"] for point in folding_points],
            mode="lines",
            name="Free energy of folding (kcal/mol)",
            line={"color": "#7c3aed", "width": 2.2},
        )
    )
    folding_figure.update_layout(
        margin={"l": 24, "r": 24, "t": 24, "b": 24},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        hovermode="x unified",
        showlegend=False,
    )
    folding_figure.update_xaxes(title_text="pH")
    folding_figure.update_yaxes(title_text="Free energy of folding (kcal/mol)")
    st.plotly_chart(
        folding_figure,
        width="stretch",
        config={"displaylogo": False},
    )
