from __future__ import annotations

from mdprotonation.viewer import compute_residue_focus_targets


def test_compute_residue_focus_targets_groups_atoms_by_residue() -> None:
    pdb_text = (
        "ATOM      1  N   GLU A  12      10.000  11.000  12.000  1.00 20.00           N\n"
        "ATOM      2  CA  GLU A  12      12.000  11.000  12.000  1.00 20.00           C\n"
        "ATOM      3  C   GLU A  12      11.000  13.000  12.000  1.00 20.00           C\n"
        "ATOM      4  N   ASP B 105      -1.000   0.000   0.000  1.00 20.00           N\n"
    )

    focus_targets = compute_residue_focus_targets(pdb_text)

    assert set(focus_targets) == {("A", 12, ""), ("B", 105, "")}
    glu_target = focus_targets[("A", 12, "")]
    assert glu_target.center == (11.0, 11.666666666666666, 12.0)
    assert glu_target.radius > 2.0
    assert glu_target.token == "A:12:-"


def test_compute_residue_focus_targets_preserves_insertion_code() -> None:
    pdb_text = (
        "ATOM      1  N   HIS A  42A      0.000   0.000   0.000  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A  42A      0.000   0.000   2.000  1.00 20.00           C\n"
    )

    focus_targets = compute_residue_focus_targets(pdb_text)

    his_target = focus_targets[("A", 42, "A")]
    assert his_target.token == "A:42:A"
    assert his_target.center == (0.0, 0.0, 1.0)
