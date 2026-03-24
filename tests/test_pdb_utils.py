from __future__ import annotations

from pkaScope.pdb_utils import parse_pdb_atoms


def test_parse_pdb_atoms_returns_residue_keys_and_coordinates() -> None:
    pdb_text = (
        "ATOM      1  N   HIS A  42A      10.000  11.000  12.500  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A  42A      10.500  11.500  13.000  1.00 20.00           C\n"
    )

    atoms = parse_pdb_atoms(pdb_text, filename="test.pdb")

    assert len(atoms) == 2
    assert atoms[0].residue_key == ("A", 42, "A")
    assert atoms[0].coordinates == (10.0, 11.0, 12.5)
