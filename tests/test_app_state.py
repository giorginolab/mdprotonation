from __future__ import annotations

from pkaScope.app_state import ph7_reference_state


def test_ph7_reference_state_uses_requested_residue_classes() -> None:
    assert ph7_reference_state("ARG") == "Cationic"
    assert ph7_reference_state("HIS") == "Cationic"
    assert ph7_reference_state("LYS") == "Cationic"
    assert ph7_reference_state("ASP") == "Anionic"
    assert ph7_reference_state("GLU") == "Anionic"
    assert ph7_reference_state("CYS") == "Neutral"


def test_ph7_reference_state_normalizes_case_and_spacing() -> None:
    assert ph7_reference_state("  his ") == "Cationic"
    assert ph7_reference_state(" glu") == "Anionic"
