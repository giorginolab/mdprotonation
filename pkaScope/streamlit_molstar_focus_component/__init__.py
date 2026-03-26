from __future__ import annotations

from pathlib import Path

import streamlit.components.v1 as components

_COMPONENT_FUNC = components.declare_component(
    "pkaScope_molstar_focus_component",
    path=str(Path(__file__).resolve().parent / "frontend"),
)


def st_molstar_focus_component(
    *,
    mvs_data: str,
    model_key: str,
    focus_residue: dict[str, object] | None,
    height: int | None,
    focus_duration_ms: int = 420,
    key: str | None = None,
):
    return _COMPONENT_FUNC(
        mvsData=mvs_data,
        mvsFormat="mvsx",
        modelKey=model_key,
        focusResidue=focus_residue,
        focusDurationMs=focus_duration_ms,
        height=height,
        key=key,
        default=None,
    )
