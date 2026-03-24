from __future__ import annotations

import streamlit as st

from ..app_state import AppState
from ..pka_plot import create_pka_comparison_plot_figure, create_pka_plot_figure
from ..presentation import PKA_PLOT_CAPTION


def render_pka_plot_tab(app_state: AppState) -> None:
    st.subheader("Residue pKa Landscape")
    site_states = list(app_state.site_states)
    pka_plot_figure = create_pka_plot_figure(site_states, app_state.ph)
    st.plotly_chart(
        pka_plot_figure,
        use_container_width=True,
        config={"displaylogo": False},
    )
    st.caption(PKA_PLOT_CAPTION)

    st.divider()
    st.subheader("Computed vs Model pKa")
    comparison_figure = create_pka_comparison_plot_figure(site_states, app_state.ph)
    st.pyplot(comparison_figure, use_container_width=True)
    st.caption("Comparison between computed pKa (blue) and standard model pKa (orange).")
