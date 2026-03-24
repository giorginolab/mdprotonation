from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

SelectionSource = Literal["main", "responsive"]


@dataclass(frozen=True)
class TableSelection:
    rows: tuple[int, ...]

    @property
    def first_row(self) -> int | None:
        if not self.rows:
            return None
        return self.rows[0]


def parse_selected_rows(selection: object | None) -> tuple[int, ...]:
    if selection is None:
        return ()

    rows = getattr(selection, "rows", None)
    if rows is None and isinstance(selection, dict):
        rows = selection.get("rows")
    if isinstance(rows, list):
        parsed_rows = tuple(row for row in rows if isinstance(row, int))
        if parsed_rows:
            return parsed_rows

    cells = getattr(selection, "cells", None)
    if cells is None and isinstance(selection, dict):
        cells = selection.get("cells")
    if not isinstance(cells, list):
        return ()

    unique_rows: list[int] = []
    for cell in cells:
        if isinstance(cell, (tuple, list)) and cell and isinstance(cell[0], int):
            if cell[0] not in unique_rows:
                unique_rows.append(cell[0])
    return tuple(unique_rows)


def get_table_selection(key: str, event_selection: object | None = None) -> TableSelection:
    import streamlit as st

    event_rows = parse_selected_rows(event_selection)
    if event_rows:
        return TableSelection(event_rows)

    state = st.session_state.get(key)
    if state is None:
        return TableSelection(())

    selection = getattr(state, "selection", None)
    if selection is None and isinstance(state, dict):
        selection = state.get("selection")
    return TableSelection(parse_selected_rows(selection))


def set_active_selection_source(source: SelectionSource) -> None:
    import streamlit as st

    st.session_state["active_selection_source"] = source


def get_active_selection_source(default: SelectionSource = "main") -> SelectionSource:
    import streamlit as st

    source = st.session_state.get("active_selection_source", default)
    if source in ("main", "responsive"):
        return source
    return default
