from __future__ import annotations

import streamlit as st

from ..app_state import AppState
from ..pka_plot import create_pka_plot_figure
from ..presentation import PKA_PLOT_CAPTION


def render_pka_plot_tab(app_state: AppState) -> None:
    st.subheader("Residue pKa Landscape")
    pka_plot_figure = create_pka_plot_figure(list(app_state.site_states), app_state.ph)
    st.pyplot(pka_plot_figure, use_container_width=True)
    st.caption(PKA_PLOT_CAPTION)
