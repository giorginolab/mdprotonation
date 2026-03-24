from __future__ import annotations

from functools import partial

import pandas as pd
import streamlit as st

from ..app_state import AppState
from ..protonation import SiteState
from ..ui_state import get_table_selection, set_active_selection_source


def render_sidebar_summary(app_state: AppState) -> SiteState | None:
    with st.sidebar:
        st.subheader("Current State")
        metric_col1, metric_col2 = st.columns(2)
        metric_col1.metric("Titratable sites", len(app_state.site_states))
        metric_col2.metric("Transitioning", app_state.transitioning_count)
        metric_col1.metric("Folded charge", f"{app_state.folded_charge:.2f}")
        metric_col2.metric("Structure", app_state.analysis.structure_name)

        st.markdown("#### Residues responding most at this pH")
        responsive_df = pd.DataFrame(
            [
                {
                    "Site": state.site.label,
                    "pKa": state.site.pka,
                    "State": state.dominant_state,
                    "% Protonated": state.protonated_fraction * 100.0,
                }
                for state in app_state.top_responsive_sites
            ]
        )
        table_event = st.dataframe(
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
        selection = get_table_selection("responsive-site-table", table_event.selection)
        if selection.first_row is None:
            return None
        if selection.first_row >= len(app_state.top_responsive_sites):
            return None
        return app_state.top_responsive_sites[selection.first_row]
