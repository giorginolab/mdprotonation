from __future__ import annotations

from pkaScope.pdb_utils import load_pdb_universe, parse_pdb_atoms, preprocess_pdb_text


def test_parse_pdb_atoms_returns_residue_keys_and_coordinates() -> None:
    pdb_text = (
        "ATOM      1  N   HIS A  42A     10.000  11.000  12.500  1.00 20.00           N\n"
        "ATOM      2  CA  HIS A  42A     10.500  11.500  13.000  1.00 20.00           C\n"
    )

    atoms = parse_pdb_atoms(pdb_text, filename="test.pdb")

    assert len(atoms) == 2
    assert atoms[0].residue_key == ("A", 42, "A")
    assert atoms[0].coordinates == (10.0, 11.0, 12.5)


def test_preprocess_pdb_text_renames_histidine_variants_to_his() -> None:
    pdb_text = (
        "ATOM      1  N   HID A  10      10.000  11.000  12.500  1.00 20.00           N\n"
        "ATOM      2  CA AHIE A  11      10.500  11.500  13.000  0.60 20.00           C\n"
        "ATOM      3  CA BHIE A  11      10.700  11.700  13.200  0.40 20.00           C\n"
        "ATOM      4  H   HSP A  12      11.000  12.000  13.500  1.00 20.00           H\n"
        "ATOM      5  C   HSP A  12      11.200  12.200  13.700  1.00 20.00           C\n"
    )

    processed = preprocess_pdb_text(pdb_text)
    universe = load_pdb_universe(processed)

    assert " HID " not in processed
    assert " HIE " not in processed
    assert " HSP " not in processed
    assert " HIS " in processed
    assert all(not atom.name.strip().upper().startswith("H") for atom in universe.atoms)
    ca_atoms = universe.select_atoms("resid 11 and name CA")
    assert len(ca_atoms) == 1
