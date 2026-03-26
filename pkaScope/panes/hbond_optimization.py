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
HBOND_SELECTED_RESIDUE_KEY = "hbond-diagnostics-selected-residue"

_EXPLICIT_RESIDUE_PATTERNS = (
    re.compile(r"\bin residue\s+([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)\b"),
    re.compile(r"\bUnable to debump\s+([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)\b"),
)
_RESNUM_CHAIN_PATTERN = re.compile(
    r"^\s*([A-Z0-9]{2,3})\s+(-?\d+)([A-Za-z]?)\s+([A-Za-z0-9])\s*$"
)
_GENERAL_RESIDUE_PATTERN = re.compile(
    r"\b([A-Z0-9]{2,3})\s+([A-Za-z0-9])\s+(-?\d+)([A-Za-z]?)\b"
)
_UNKNOWN_RESIDUE_HEADER_PATTERN = re.compile(
    r"could not identify the following residues and residue numbers",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class DiagnosticResidueMention:
    line_index: int
    line_text: str
    category: str
    residue_type: str
    chain_id: str
    residue_number: int
    insertion_code: str

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain_id, self.residue_number, self.insertion_code)

    @property
    def residue_label(self) -> str:
        residue_token = (
            f"{self.residue_number:>5}{self.insertion_code}"
            if self.insertion_code
            else f"{self.residue_number:>5}"
        )
        return f"{self.chain_id}: {residue_token} {self.residue_type}"


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
        "Jurrus E, et al. Improvements to the APBS biomolecular solvation "
        "software suite. Protein Sci 27 112-128 (2018)."
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
    diagnostics_lines = diagnostics_preview_lines(result.diagnostics)
    mentions = extract_diagnostic_residue_mentions(diagnostics_lines)
    selected_focus_key = _selected_focus_key_from_state()
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

    _render_clickable_diagnostics(diagnostics_lines, mentions)

    if selected_focus_key is not None and before_focus is None and after_focus is None:
        st.caption(
            "Selected residue from diagnostics could not be located in the current structures."
        )


def diagnostics_preview_lines(
    diagnostics: tuple[str, ...],
) -> tuple[str, ...]:
    return diagnostics


def _render_clickable_diagnostics(
    diagnostics_lines: tuple[str, ...],
    mentions: tuple[DiagnosticResidueMention, ...],
) -> None:
    st.markdown("#### Diagnostics")
    if not diagnostics_lines:
        st.caption("No diagnostics were produced by pdb2pqr.")
        return

    sorted_mentions = sort_mentions_natural(mentions)
    if sorted_mentions:
        st.caption("Click a residue mention to focus both viewers on that residue.")
        mention_df = pd.DataFrame(
            [
                {
                    "Category": mention.category,
                    "Residue": mention.residue_label,
                    "Message": mention.line_text,
                }
                for mention in sorted_mentions
            ]
        )
        table_event = st.dataframe(
            mention_df,
            width="stretch",
            hide_index=True,
            key="hbond-diagnostics-table",
            selection_mode="single-row",
            on_select="rerun",
        )
        selection = get_table_selection("hbond-diagnostics-table", table_event.selection)
        selected_row = selection.first_row
        if selected_row is not None and selected_row < len(sorted_mentions):
            selected_residue_key = sorted_mentions[selected_row].residue_key
            if st.session_state.get(HBOND_SELECTED_RESIDUE_KEY) != selected_residue_key:
                st.session_state[HBOND_SELECTED_RESIDUE_KEY] = selected_residue_key
                st.rerun()
            st.caption(
                f"Selected diagnostic residue: `{sorted_mentions[selected_row].residue_label}`"
            )
        elif HBOND_SELECTED_RESIDUE_KEY in st.session_state:
            del st.session_state[HBOND_SELECTED_RESIDUE_KEY]
    else:
        st.caption("No residue mentions detected in diagnostics.")

    with st.expander("Full diagnostics log", expanded=False):
        st.code("\n".join(diagnostics_lines), language="text")


def _selected_focus_key_from_state() -> tuple[str, int, str] | None:
    selected_focus_key = st.session_state.get(HBOND_SELECTED_RESIDUE_KEY)
    if not isinstance(selected_focus_key, tuple) or len(selected_focus_key) != 3:
        return None
    chain_id, residue_number, insertion_code = selected_focus_key
    if not isinstance(chain_id, str) or not isinstance(residue_number, int):
        return None
    if not isinstance(insertion_code, str):
        return None
    return (chain_id, residue_number, insertion_code)


def extract_diagnostic_residue_mentions(
    diagnostics: tuple[str, ...],
) -> tuple[DiagnosticResidueMention, ...]:
    mentions: list[DiagnosticResidueMention] = []
    in_unknown_residue_block = False
    for line_index, line_text in enumerate(diagnostics):
        if _UNKNOWN_RESIDUE_HEADER_PATTERN.search(line_text):
            in_unknown_residue_block = True
            continue
        category = categorize_diagnostic_line(line_text)
        if in_unknown_residue_block and category == "Other":
            category = "Unknown residue"
        seen: set[tuple[int, int, int, int]] = set()
        patterns = list(_EXPLICIT_RESIDUE_PATTERNS)
        if in_unknown_residue_block:
            patterns.append(_RESNUM_CHAIN_PATTERN)
        if not any(pattern.search(line_text) for pattern in _EXPLICIT_RESIDUE_PATTERNS):
            patterns.append(_GENERAL_RESIDUE_PATTERN)
        matched_line = False
        for pattern in patterns:
            for match in pattern.finditer(line_text):
                match_key = (match.start(), match.end(), line_index, hash(pattern.pattern))
                if match_key in seen:
                    continue
                seen.add(match_key)
                residue_type = match.group(1).strip().upper()
                if pattern is _RESNUM_CHAIN_PATTERN:
                    try:
                        residue_number = int(match.group(2))
                    except ValueError:
                        continue
                    insertion_code = match.group(3).strip().upper()
                    chain_id = match.group(4).strip() or "?"
                else:
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
                        category=category,
                        residue_type=residue_type,
                        chain_id=chain_id,
                        residue_number=residue_number,
                        insertion_code=insertion_code,
                    )
                )
                matched_line = True
        if in_unknown_residue_block and not matched_line:
            in_unknown_residue_block = False
    return tuple(mentions)


def categorize_diagnostic_line(line_text: str) -> str:
    lowered = line_text.lower()
    if "missing atom" in lowered:
        return "Missing atom"
    if "unable to debump" in lowered:
        return "Debump"
    if "water optimization" in lowered:
        return "Water optimization"
    if "unable to find amino or nucleic acid definition" in lowered:
        return "Unknown residue"
    return "Other"


def sort_mentions_natural(
    mentions: tuple[DiagnosticResidueMention, ...],
) -> tuple[DiagnosticResidueMention, ...]:
    return tuple(
        sorted(
            mentions,
            key=lambda mention: (
                _chain_natural_sort_key(mention.chain_id),
                mention.residue_number,
                mention.insertion_code,
                mention.residue_type,
                mention.line_index,
            ),
        )
    )


def _chain_natural_sort_key(chain_id: str) -> tuple[int, int | str]:
    if chain_id.isdigit():
        return (0, int(chain_id))
    if len(chain_id) == 1 and chain_id.isalpha():
        return (1, ord(chain_id.upper()))
    return (2, chain_id.upper())
