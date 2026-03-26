from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import streamlit as st

from ..pdb2pqr_analysis import (
    Pdb2pqrAnalysisError,
    Pdb2pqrRunOptions,
    Pdb2pqrRunResult,
    run_pdb2pqr_analysis,
)
from ..propka_analysis import PropkaAnalysis
from ..ui_state import get_table_selection
from ..viewer import compute_residue_focus_targets, st_molstar_focusable_content

HBOND_SCHEMA_VERSION = "2026-03-26-hbond-opt-v1"
HBOND_LAST_REQUEST_KEY = "hbond-last-request"
HBOND_LAST_RESULT_KEY = "hbond-last-result"
HBOND_VIEWER_HEIGHT = 520

_EXPLICIT_RESIDUE_PATTERNS = (
    re.compile(r"\bin residue\s+([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)\b"),
    re.compile(r"\bUnable to debump\s+([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)\b"),
)
_GENERAL_RESIDUE_PATTERN = re.compile(
    r"\b([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)\b"
)


@dataclass(frozen=True)
class DiagnosticResidueMention:
    line_index: int
    line_text: str
    residue_type: str
    chain_id: str
    residue_number: int
    insertion_code: str

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain_id, self.residue_number, self.insertion_code)

    @property
    def residue_label(self) -> str:
        insertion = self.insertion_code if self.insertion_code else ""
        return f"{self.residue_type} {self.chain_id} {self.residue_number}{insertion}"


@st.cache_data(show_spinner=False)
def cached_hbond_optimization(
    pdb_text: str,
    source_name: str,
    ph: float,
    schema_version: str,
) -> Pdb2pqrRunResult:
    del schema_version
    return run_pdb2pqr_analysis(
        pdb_text,
        source_name,
        Pdb2pqrRunOptions(ph=ph),
    )


@st.cache_data(show_spinner=False)
def cached_focus_targets(pdb_text: str) -> dict[tuple[str, int, str], object]:
    return compute_residue_focus_targets(pdb_text)


def render_hbond_optimization_tab(analysis: PropkaAnalysis, ph: float) -> None:
    st.subheader("H-Bond Placement and Optimization (pdb2pqr)")
    st.caption(
        "Run pdb2pqr to place hydrogens and optimize hydrogen-bond orientations."
    )
    st.markdown("#### Citations")
    st.caption(
        "Jurrus E, et al. Improvements to the APBS biomolecular solvation "
        "software suite. Protein Sci 27 112-128 (2018)."
    )
    st.caption(
        "Dolinsky TJ, et al. PDB2PQR: expanding and upgrading automated "
        "preparation of biomolecular structures for molecular simulations. "
        "Nucleic Acids Res 35 W522-W525 (2007)."
    )

    options = Pdb2pqrRunOptions(ph=ph)
    current_signature = request_signature(
        source_name=analysis.source_name,
        ph=ph,
    )
    _render_settings(options)

    run_clicked = st.button(
        "Run pdb2pqr optimization",
        type="primary",
        key="hbond-opt-run-button",
    )
    if run_clicked:
        try:
            with st.spinner("Running pdb2pqr..."):
                result = cached_hbond_optimization(
                    analysis.pdb_text,
                    analysis.source_name,
                    ph,
                    HBOND_SCHEMA_VERSION,
                )
        except Pdb2pqrAnalysisError as exc:
            st.error(exc.user_message)
            if exc.details:
                with st.expander("Error details"):
                    st.code(exc.details, language="text")
            return
        except Exception as exc:
            st.error("Unexpected error while running pdb2pqr.")
            st.exception(exc)
            return
        st.session_state[HBOND_LAST_REQUEST_KEY] = current_signature
        st.session_state[HBOND_LAST_RESULT_KEY] = result

    last_signature = st.session_state.get(HBOND_LAST_REQUEST_KEY)
    last_result = st.session_state.get(HBOND_LAST_RESULT_KEY)

    if last_signature is None or last_result is None:
        st.info("Click `Run pdb2pqr optimization` to generate optimized outputs.")
        return

    if not signatures_match(last_signature, current_signature):
        st.info(
            "Input structure or pH changed. Click `Run pdb2pqr optimization` "
            "to refresh results."
        )
        return

    render_result_blocks(analysis, last_result)


def _render_settings(options: Pdb2pqrRunOptions) -> None:
    st.markdown("#### Current Settings")
    settings_col1, settings_col2, settings_col3 = st.columns(3)
    settings_col1.metric("pH (global slider)", f"{options.ph:.2f}")
    settings_col2.metric("Force field", options.force_field)
    settings_col3.metric("H-bond optimization", "Enabled" if options.optimize_hbonds else "Disabled")
    st.caption("Defaults: keep-chain enabled, titration state method = propka.")


def request_signature(*, source_name: str, ph: float) -> tuple[str, float]:
    return (source_name, round(ph, 3))


def signatures_match(
    previous_signature: tuple[str, float],
    current_signature: tuple[str, float],
) -> bool:
    return previous_signature == current_signature


def render_result_blocks(analysis: PropkaAnalysis, result: Pdb2pqrRunResult) -> None:
    before_focus_targets = cached_focus_targets(analysis.pdb_text)
    after_focus_targets = cached_focus_targets(result.optimized_pdb_text)
    selected_focus_key = _render_clickable_diagnostics(result.diagnostics)
    before_focus = (
        before_focus_targets.get(selected_focus_key)
        if selected_focus_key is not None
        else None
    )
    after_focus = (
        after_focus_targets.get(selected_focus_key)
        if selected_focus_key is not None
        else None
    )

    st.markdown("#### Run Summary")
    metric_col1, metric_col2, metric_col3 = st.columns(3)
    metric_col1.metric("Atoms in optimized structure", f"{result.atom_count}")
    metric_col2.metric("Residues in optimized structure", f"{result.residue_count}")
    metric_col3.metric("Warnings in diagnostics", f"{result.warning_count}")

    st.markdown("#### Structure Preview")
    before_col, after_col = st.columns(2, gap="large")
    with before_col:
        st.caption("Before optimization")
        st_molstar_focusable_content(
            analysis.pdb_text,
            "pdb",
            file_name=f"{Path(analysis.source_name).stem}-original.pdb",
            focus_residue=before_focus,
            height=f"{HBOND_VIEWER_HEIGHT}px",
            key=f"hbond-before-{analysis.structure_name}-{analysis.source_name}",
        )
    with after_col:
        st.caption("After optimization")
        st_molstar_focusable_content(
            result.optimized_pdb_text,
            "pdb",
            file_name=f"{Path(analysis.source_name).stem}-optimized.pdb",
            focus_residue=after_focus,
            height=f"{HBOND_VIEWER_HEIGHT}px",
            key=f"hbond-after-{analysis.structure_name}-{analysis.source_name}",
        )

    st.markdown("#### Download Outputs")
    download_col1, download_col2 = st.columns(2)
    with download_col1:
        st.download_button(
            "Download optimized PDB",
            data=result.optimized_pdb_text,
            file_name=f"{Path(analysis.source_name).stem}-optimized.pdb",
            mime="chemical/x-pdb",
            key=f"hbond-download-pdb-{analysis.source_name}",
        )
    with download_col2:
        if result.pqr_text is None:
            st.caption("PQR output unavailable for this run.")
        else:
            st.download_button(
                "Download PQR",
                data=result.pqr_text,
                file_name=f"{Path(analysis.source_name).stem}.pqr",
                mime="text/plain",
                key=f"hbond-download-pqr-{analysis.source_name}",
            )

    if selected_focus_key is not None and before_focus is None and after_focus is None:
        st.caption(
            "Selected residue from diagnostics could not be located in the current structures."
        )


def diagnostics_preview_lines(
    diagnostics: tuple[str, ...],
) -> tuple[str, ...]:
    return diagnostics


def _render_clickable_diagnostics(
    diagnostics: tuple[str, ...],
) -> tuple[str, int, str] | None:
    st.markdown("#### Diagnostics")
    diagnostics_lines = diagnostics_preview_lines(diagnostics)
    if not diagnostics_lines:
        st.caption("No diagnostics were produced by pdb2pqr.")
        return None

    mentions = extract_diagnostic_residue_mentions(diagnostics_lines)
    selected_focus_key: tuple[str, int, str] | None = None
    if mentions:
        st.caption("Click a residue mention to focus both viewers on that residue.")
        mention_df = pd.DataFrame(
            [
                {
                    "Residue": mention.residue_label,
                    "Line": mention.line_index + 1,
                    "Message": mention.line_text,
                }
                for mention in mentions
            ]
        )
        table_event = st.dataframe(
            mention_df,
            width="stretch",
            hide_index=True,
            key="hbond-diagnostics-table",
            selection_mode="single-row",
        )
        selection = get_table_selection("hbond-diagnostics-table", table_event.selection)
        selected_row = selection.first_row
        if selected_row is not None and selected_row < len(mentions):
            selected_focus_key = mentions[selected_row].residue_key
            st.caption(f"Selected diagnostic residue: `{mentions[selected_row].residue_label}`")
    else:
        st.caption("No residue mentions detected in diagnostics.")

    with st.expander("Full diagnostics log", expanded=False):
        st.code("\n".join(diagnostics_lines), language="text")
    return selected_focus_key


def extract_diagnostic_residue_mentions(
    diagnostics: tuple[str, ...],
) -> tuple[DiagnosticResidueMention, ...]:
    mentions: list[DiagnosticResidueMention] = []
    for line_index, line_text in enumerate(diagnostics):
        seen: set[tuple[int, int, int, int]] = set()
        patterns = list(_EXPLICIT_RESIDUE_PATTERNS)
        if not any(pattern.search(line_text) for pattern in _EXPLICIT_RESIDUE_PATTERNS):
            patterns.append(_GENERAL_RESIDUE_PATTERN)
        for pattern in patterns:
            for match in pattern.finditer(line_text):
                match_key = (match.start(), match.end(), line_index, hash(pattern.pattern))
                if match_key in seen:
                    continue
                seen.add(match_key)
                residue_type = match.group(1).strip().upper()
                chain_id = match.group(2).strip() or "?"
                try:
                    residue_number = int(match.group(3))
                except ValueError:
                    continue
                insertion_code = match.group(4).strip().upper()
                mentions.append(
                    DiagnosticResidueMention(
                        line_index=line_index,
                        line_text=line_text,
                        residue_type=residue_type,
                        chain_id=chain_id,
                        residue_number=residue_number,
                        insertion_code=insertion_code,
                    )
                )
    return tuple(mentions)
