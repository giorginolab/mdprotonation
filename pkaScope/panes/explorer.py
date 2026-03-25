from __future__ import annotations

from functools import partial

import pandas as pd
import streamlit as st

from ..app_state import AppState
from ..presentation import (
    CHARGE_SCALE_CAPTION,
    PH7_STATE_CAPTION,
    TRANSITION_MARKER,
    VIEWER_ENCODING_CAPTION,
    style_charge_row,
)
from ..protonation import SiteState
from ..ui_state import (
    get_active_selection_source,
    get_table_selection,
    set_active_selection_source,
)
from ..viewer import st_molstar_focusable_content

EXPLORER_VIEWER_HEIGHT = 640
SITE_TABLE_MIN_HEIGHT = 380
SITE_TABLE_MAX_HEIGHT = 640
SITE_TABLE_HEADER_HEIGHT = 38
SITE_TABLE_ROW_HEIGHT = 35
SITE_TABLE_VISIBLE_COLUMNS = (
    "pH7",
    "Site",
    "pKa",
    "Model pKa",
    "% Buried",
    "% Prot",
    "Charge",
    "State",
)


def render_explorer_tab(
    app_state: AppState,
    transition_only: bool,
    responsive_selected_state: SiteState | None,
) -> None:
    sorted_states = app_state.filtered_site_states(transition_only)
    main_selected_state = _selected_state_from_index(
        sorted_states,
        get_table_selection("site-state-table").first_row,
    )
    selected_state = _resolve_selected_state(
        main_selected_state=main_selected_state,
        responsive_selected_state=responsive_selected_state,
    )
    focus_residue = (
        app_state.residue_focus_targets.get(selected_state.site.residue_key)
        if selected_state is not None
        else None
    )
    state_df = _build_site_state_dataframe(app_state, sorted_states)
    styled_state_df = state_df.style.apply(style_charge_row, axis=1).format(
        {
            "pKa": "{:.2f}",
            "Model pKa": "{:.2f}",
            "Buried fraction": "{:.3f}",
            "% Buried": "{:.0f}",
            "Charge": "{:.3f}",
            "% Prot": "{:.0f}",
        }
    )

    st.caption(
        "Split view keeps the residue list and the 3D structure visible together. "
        "Select a residue row to refocus the camera."
    )
    viewer_col, table_col = st.columns((5, 4), gap="large")

    with viewer_col:
        st.subheader("Structure Viewer")
        st_molstar_focusable_content(
            app_state.encoded_pdb_text,
            "pdb",
            file_name=f"{app_state.analysis.structure_name}-ph-{app_state.ph:.1f}.pdb",
            focus_residue=focus_residue,
            height=f"{EXPLORER_VIEWER_HEIGHT}px",
            key=f"molstar-{app_state.analysis.structure_name}-{app_state.ph:.1f}",
        )
        st.caption(VIEWER_ENCODING_CAPTION)
        if selected_state is not None:
            st.caption(
                f"Camera focus target: `{selected_state.site.label}` from the selected table row."
            )

    with table_col:
        st.subheader("Residue Protonation States")
        table_event = st.dataframe(
            styled_state_df,
            width="stretch",
            hide_index=True,
            height=_site_table_height(len(state_df)),
            column_order=SITE_TABLE_VISIBLE_COLUMNS,
            key="site-state-table",
            on_select=partial(set_active_selection_source, "main"),
            selection_mode="single-cell",
        )
        st.caption(CHARGE_SCALE_CAPTION)
        st.caption(PH7_STATE_CAPTION)

    table_selection = get_table_selection("site-state-table", table_event.selection)
    main_selected_state = _selected_state_from_index(sorted_states, table_selection.first_row)
    if main_selected_state is not None and get_active_selection_source() == "main":
        selected_state = main_selected_state
    elif selected_state is None:
        selected_state = responsive_selected_state

    st.subheader("Selected Site Interactions")
    if selected_state is None:
        st.info("Click a row in either residue table to inspect PROPKA interaction residues.")
        return

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
        st.dataframe(interaction_rows, width="stretch", hide_index=True)
    else:
        st.info("PROPKA reported no non-self interaction residues for this site.")


def _build_site_state_dataframe(
    app_state: AppState,
    sorted_states: tuple[SiteState, ...],
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "pH7": (
                    ""
                    if state.dominant_state
                    == app_state.standard_state_by_site.get(
                        state.site.site_key,
                        state.dominant_state,
                    )
                    else TRANSITION_MARKER
                ),
                "Site": _site_column_label(state),
                "Residue": state.site.residue_type,
                "Chain": state.site.chain_id,
                "Residue #": state.site.residue_number,
                "pKa": state.site.pka,
                "Model pKa": state.site.model_pka,
                "Buried fraction": state.site.buried_fraction,
                "% Buried": state.site.buried_fraction * 100.0,
                "% Prot": state.protonated_fraction * 100.0,
                "Charge": state.current_charge,
                "State": state.dominant_state,
            }
            for state in sorted_states
        ]
    )


def _resolve_selected_state(
    *,
    main_selected_state: SiteState | None,
    responsive_selected_state: SiteState | None,
) -> SiteState | None:
    if (
        get_active_selection_source() == "responsive"
        and responsive_selected_state is not None
    ):
        return responsive_selected_state
    if main_selected_state is not None:
        return main_selected_state
    return responsive_selected_state


def _selected_state_from_index(
    states: tuple[SiteState, ...],
    row_index: int | None,
) -> SiteState | None:
    if row_index is None or row_index >= len(states):
        return None
    return states[row_index]


def _site_table_height(row_count: int) -> int:
    if row_count <= 0:
        return SITE_TABLE_MIN_HEIGHT
    estimated_height = SITE_TABLE_HEADER_HEIGHT + (row_count * SITE_TABLE_ROW_HEIGHT)
    return max(SITE_TABLE_MIN_HEIGHT, min(SITE_TABLE_MAX_HEIGHT, estimated_height))


def _site_column_label(state: SiteState) -> str:
    insertion_code = state.site.insertion_code.strip()
    residue_token = (
        f"{state.site.residue_number:>5}{insertion_code}"
        if insertion_code
        else f"{state.site.residue_number:>5}"
    )
    return f"{state.site.chain_id}:{residue_token}-{state.site.residue_type}"
