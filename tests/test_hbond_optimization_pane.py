from __future__ import annotations

from pkaScope.panes.hbond_optimization import (
    diagnostics_preview_lines,
    extract_diagnostic_residue_mentions,
    request_signature,
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
