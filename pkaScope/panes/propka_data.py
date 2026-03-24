from __future__ import annotations

import streamlit as st

from ..propka_analysis import PropkaAnalysis


def render_propka_data_tab(analysis: PropkaAnalysis) -> None:
    st.subheader("PROPKA Summary")
    st.code(analysis.summary_text.strip(), language="text")
    st.subheader("PROPKA Determinants")
    st.code(analysis.determinants_text.strip(), language="text")
