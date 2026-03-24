from __future__ import annotations

from mdprotonation.ui_state import parse_selected_rows


def test_parse_selected_rows_prefers_explicit_row_selection() -> None:
    selection = {"rows": [3, 1, "x"], "cells": [[4, "Site"]]}

    assert parse_selected_rows(selection) == (3, 1)


def test_parse_selected_rows_deduplicates_rows_from_selected_cells() -> None:
    selection = {"cells": [[2, "Site"], [2, "pKa"], [4, "State"]]}

    assert parse_selected_rows(selection) == (2, 4)
