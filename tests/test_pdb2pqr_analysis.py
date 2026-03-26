from __future__ import annotations

import pytest

from pkaScope.pdb2pqr_analysis import (
    Pdb2pqrAnalysisError,
    Pdb2pqrRunOptions,
    run_pdb2pqr_analysis,
)


def test_run_pdb2pqr_analysis_maps_success_payload(monkeypatch: pytest.MonkeyPatch) -> None:
    def _fake_execute_pdb2pqr_python_api(**_: object) -> dict[str, object]:
        return {
            "optimized_pdb_text": (
                "ATOM      1  N   GLY A   1      11.100  11.200  11.300  1.00 10.00           N\n"
                "ATOM      2  CA  GLY A   1      12.100  11.200  11.300  1.00 10.00           C\n"
                "ATOM      3  N   ASP A   2      13.100  11.200  11.300  1.00 10.00           N\n"
            ),
            "pqr_text": "ATOM ...\n",
            "diagnostics_text": "WARNING: test warning\nAll done\n",
        }

    monkeypatch.setattr(
        "pkaScope.pdb2pqr_analysis._execute_pdb2pqr_python_api",
        _fake_execute_pdb2pqr_python_api,
    )

    result = run_pdb2pqr_analysis(
        "ATOM ...\n",
        "demo.pdb",
        Pdb2pqrRunOptions(ph=7.0),
    )

    assert result.pqr_text == "ATOM ...\n"
    assert result.atom_count == 3
    assert result.residue_count == 2
    assert result.warning_count == 1
    assert result.diagnostics == ("WARNING: test warning", "All done")


def test_run_pdb2pqr_analysis_reports_missing_dependency(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def _raise_missing_dependency(**_: object) -> dict[str, object]:
        raise ModuleNotFoundError("No module named 'pdb2pqr'")

    monkeypatch.setattr(
        "pkaScope.pdb2pqr_analysis._execute_pdb2pqr_python_api",
        _raise_missing_dependency,
    )

    with pytest.raises(Pdb2pqrAnalysisError) as exc_info:
        run_pdb2pqr_analysis(
            "ATOM ...\n",
            "demo.pdb",
            Pdb2pqrRunOptions(ph=7.0),
        )

    assert exc_info.value.code == "missing_dependency"
    assert "pdb2pqr is not installed" in exc_info.value.user_message
