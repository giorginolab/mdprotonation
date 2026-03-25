from __future__ import annotations

import pytest

from pkaScope.chain_shift import (
    compare_chain_pkas,
    compare_site_sets_pkas,
    extract_chain_pdb,
    extract_chains_pdb,
    list_chain_ids,
    top_shift_residue_keys,
)
from pkaScope.propka_analysis import TitrationSite


def _make_site(
    *,
    label: str,
    residue_type: str,
    chain_id: str,
    residue_number: int,
    pka: float,
    insertion_code: str = "",
) -> TitrationSite:
    return TitrationSite(
        label=label,
        residue_type=residue_type,
        chain_id=chain_id,
        residue_number=residue_number,
        insertion_code=insertion_code,
        pka=pka,
        model_pka=pka,
        charged_state_charge=1.0 if residue_type in {"ARG", "LYS", "HIS"} else -1.0,
        buried_fraction=0.0,
        interactions=(),
    )


def test_list_chain_ids_preserves_order_and_normalizes_blank_chain() -> None:
    pdb_text = (
        "ATOM      1  N   ASP B  12      10.000  11.000  12.000  1.00 20.00           N\n"
        "ATOM      2  N   GLU A  24      10.000  11.000  12.000  1.00 20.00           N\n"
        "ATOM      3  N   GLY    3       9.000   9.000   9.000  1.00 20.00           N\n"
        "ATOM      4  N   ASP A  30      10.000  11.000  12.000  1.00 20.00           N\n"
    )

    assert list_chain_ids(pdb_text) == ("B", "A", "?")


def test_extract_chain_pdb_filters_coordinate_records_for_selected_chain() -> None:
    pdb_text = (
        "HEADER    TEST\n"
        "ATOM      1  N   ASP A  12      10.000  11.000  12.000  1.00 20.00           N\n"
        "ATOM      2  N   GLU B  42      13.000  14.000  15.000  1.00 20.00           N\n"
        "TER\n"
        "END\n"
    )

    extracted_pdb = extract_chain_pdb(pdb_text, "B")

    assert "ATOM" in extracted_pdb
    assert "GLU B  42" in extracted_pdb
    assert "ASP A  12" not in extracted_pdb


def test_extract_chains_pdb_filters_coordinate_records_for_selected_chains() -> None:
    pdb_text = (
        "HEADER    TEST\n"
        "ATOM      1  N   ASP A  12      10.000  11.000  12.000  1.00 20.00           N\n"
        "ATOM      2  N   GLU B  42      13.000  14.000  15.000  1.00 20.00           N\n"
        "ATOM      3  N   LYS C  55      11.000  12.000  13.000  1.00 20.00           N\n"
        "TER\n"
        "END\n"
    )

    extracted_pdb = extract_chains_pdb(pdb_text, ("A", "C"))

    assert "ASP A  12" in extracted_pdb
    assert "LYS C  55" in extracted_pdb
    assert "GLU B  42" not in extracted_pdb


def test_compare_chain_pkas_matches_by_residue_key_when_chain_ids_change() -> None:
    complex_sites = (
        _make_site(
            label="ASP  10 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=10,
            pka=4.0,
        ),
        _make_site(
            label="LYS  20 A",
            residue_type="LYS",
            chain_id="A",
            residue_number=20,
            pka=10.2,
        ),
    )
    monomer_sites = (
        _make_site(
            label="ASP  10 ?",
            residue_type="ASP",
            chain_id="?",
            residue_number=10,
            pka=4.8,
        ),
        _make_site(
            label="LYS  20 ?",
            residue_type="LYS",
            chain_id="?",
            residue_number=20,
            pka=9.6,
        ),
    )

    comparison = compare_chain_pkas(
        chain_id="A",
        complex_sites=complex_sites,
        monomer_sites=monomer_sites,
    )

    assert len(comparison.shifts) == 2
    assert comparison.unmatched_complex_labels == ()
    assert comparison.unmatched_monomer_labels == ()
    assert comparison.shifts[0].delta_pka == pytest.approx(-0.8)
    assert comparison.shifts[1].delta_pka == pytest.approx(0.6)
    assert top_shift_residue_keys(comparison.shifts, minimum_absolute_delta=0.7) == (
        ("A", 10, ""),
    )


def test_compare_site_sets_pkas_compares_multiple_chains() -> None:
    complex_sites = (
        _make_site(
            label="ASP  10 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=10,
            pka=4.2,
        ),
        _make_site(
            label="GLU  30 B",
            residue_type="GLU",
            chain_id="B",
            residue_number=30,
            pka=5.4,
        ),
    )
    monomer_sites = (
        _make_site(
            label="ASP  10 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=10,
            pka=4.8,
        ),
        _make_site(
            label="GLU  30 B",
            residue_type="GLU",
            chain_id="B",
            residue_number=30,
            pka=4.9,
        ),
    )

    comparison = compare_site_sets_pkas(
        comparison_label="Advanced chain selection",
        complex_sites=complex_sites,
        monomer_sites=monomer_sites,
    )

    assert comparison.chain_id == "Advanced chain selection"
    assert len(comparison.shifts) == 2
    assert comparison.shifts[0].delta_pka == pytest.approx(-0.6)
    assert comparison.shifts[1].delta_pka == pytest.approx(0.5)


def test_compare_site_sets_pkas_sorts_shifts_by_chain_then_residue() -> None:
    complex_sites = (
        _make_site(
            label="GLU   5 B",
            residue_type="GLU",
            chain_id="B",
            residue_number=5,
            pka=5.3,
        ),
        _make_site(
            label="ASP  40 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=40,
            pka=4.4,
        ),
    )
    monomer_sites = (
        _make_site(
            label="GLU   5 B",
            residue_type="GLU",
            chain_id="B",
            residue_number=5,
            pka=5.0,
        ),
        _make_site(
            label="ASP  40 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=40,
            pka=4.8,
        ),
    )

    comparison = compare_site_sets_pkas(
        comparison_label="Complex vs Apo",
        complex_sites=complex_sites,
        monomer_sites=monomer_sites,
    )

    assert [shift.chain_id for shift in comparison.shifts] == ["A", "B"]
    assert [shift.residue_number for shift in comparison.shifts] == [40, 5]


def test_compare_site_sets_pkas_does_not_cross_match_different_chains() -> None:
    complex_sites = (
        _make_site(
            label="ASP  12 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=12,
            pka=4.2,
        ),
        _make_site(
            label="ASP  12 C",
            residue_type="ASP",
            chain_id="C",
            residue_number=12,
            pka=4.9,
        ),
    )
    monomer_sites = (
        _make_site(
            label="ASP  12 A",
            residue_type="ASP",
            chain_id="A",
            residue_number=12,
            pka=4.6,
        ),
    )

    comparison = compare_site_sets_pkas(
        comparison_label="Complex vs Apo",
        complex_sites=complex_sites,
        monomer_sites=monomer_sites,
    )

    assert len(comparison.shifts) == 1
    assert comparison.shifts[0].chain_id == "A"
    assert comparison.unmatched_complex_labels == ("ASP  12 C",)
