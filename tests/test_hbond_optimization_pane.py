from __future__ import annotations

from pkaScope.panes.hbond_optimization import (
    DiagnosticResidueMention,
    categorize_diagnostic_line,
    diagnostics_preview_lines,
    extract_diagnostic_residue_mentions,
    request_signature,
    sort_mentions_natural,
    signatures_match,
)


def test_request_signature_tracks_source_and_ph() -> None:
    assert request_signature(source_name="a.pdb", ph=7.0) != request_signature(
        source_name="a.pdb",
        ph=7.5,
    )
    assert request_signature(source_name="a.pdb", ph=7.0) != request_signature(
        source_name="b.pdb",
        ph=7.0,
    )


def test_signatures_match_only_for_identical_values() -> None:
    previous = request_signature(source_name="demo.pdb", ph=7.0)
    assert signatures_match(previous, request_signature(source_name="demo.pdb", ph=7.0))
    assert not signatures_match(
        previous,
        request_signature(source_name="demo.pdb", ph=8.0),
    )


def test_diagnostics_preview_lines_keeps_all_logs() -> None:
    lines = tuple(f"line-{i}" for i in range(5))
    assert diagnostics_preview_lines(lines) == lines


def test_extract_diagnostic_residue_mentions_matches_common_messages() -> None:
    diagnostics = (
        "Missing atom NH2 in residue ARG D 94",
        "Missing atom OXT in residue SER D 98",
        "WARNING: Unable to debump ARG A 308",
        "Skipped atom during water optimization: H2 in HOH A 29 skipped when optimizing H2 in HOH A 8",
        "Unable to find amino or nucleic acid definition for ZN.  Parsing as new residue.",
    )

    mentions = extract_diagnostic_residue_mentions(diagnostics)
    keys = [mention.residue_key for mention in mentions]
    labels = [mention.residue_label for mention in mentions]

    assert ("D", 94, "") in keys
    assert ("D", 98, "") in keys
    assert ("A", 308, "") in keys
    assert ("A", 29, "") in keys
    assert ("A", 8, "") in keys
    assert "ZN ? 0" not in labels


def test_extract_diagnostic_residue_mentions_parses_resnum_chain_block() -> None:
    diagnostics = (
        "PDB2PQR could not identify the following residues and residue numbers as returned by PROPKA or PDB2PKA",
        "BEN 1 A",
    )
    mentions = extract_diagnostic_residue_mentions(diagnostics)
    assert len(mentions) == 1
    assert mentions[0].residue_key == ("A", 1, "")
    assert mentions[0].category == "Unknown residue"


def test_extract_diagnostic_residue_mentions_parses_unknown_definition_line() -> None:
    diagnostics = (
        "Unable to find amino or nucleic acid definition for NAG.  Parsing as new residue.",
    )
    mentions = extract_diagnostic_residue_mentions(diagnostics)
    assert len(mentions) == 1
    assert mentions[0].residue_type == "NAG"
    assert mentions[0].residue_key == ("?", 0, "")
    assert mentions[0].category == "Unknown residue"


def test_extract_diagnostic_residue_mentions_parses_failed_protonation_line() -> None:
    diagnostics = (
        "Missing atoms or failed protonation for ASN 139 A (AMD) -- please check the structure",
    )
    mentions = extract_diagnostic_residue_mentions(diagnostics)
    assert len(mentions) == 1
    assert mentions[0].residue_key == ("A", 139, "")
    assert mentions[0].residue_type == "ASN"
    assert mentions[0].category == "Failed protonation"


def test_categorize_diagnostic_line_maps_expected_categories() -> None:
    assert categorize_diagnostic_line("Missing atom NH2 in residue ARG D 94") == "Missing atom"
    assert categorize_diagnostic_line("WARNING: Unable to debump ARG A 308") == "Debump"
    assert (
        categorize_diagnostic_line(
            "Skipped atom during water optimization: H2 in HOH A 29 skipped when optimizing H2 in HOH A 8"
        )
        == "Water optimization"
    )
    assert (
        categorize_diagnostic_line(
            "Unable to find amino or nucleic acid definition for ZN.  Parsing as new residue."
        )
        == "Unknown residue"
    )
    assert (
        categorize_diagnostic_line(
            "Missing atoms or failed protonation for ASN 139 A (AMD) -- please check the structure"
        )
        == "Failed protonation"
    )


def test_sort_mentions_natural_uses_chain_then_residue_number() -> None:
    mentions = extract_diagnostic_residue_mentions(
        (
            "Missing atom NH2 in residue ARG B 10",
            "Missing atom NH2 in residue ARG B 2",
            "Missing atom NH2 in residue ARG A 5",
        )
    )
    sorted_mentions = sort_mentions_natural(mentions)
    assert [mention.residue_key for mention in sorted_mentions] == [
        ("A", 5, ""),
        ("B", 2, ""),
        ("B", 10, ""),
    ]


def test_residue_label_uses_padded_style() -> None:
    mention = DiagnosticResidueMention(
        line_index=0,
        line_text="Missing atom NH2 in residue ARG A 305",
        category="Missing atom",
        residue_type="ARG",
        chain_id="A",
        residue_number=305,
        insertion_code="",
    )
    assert mention.residue_label == "A:   305 ARG"
