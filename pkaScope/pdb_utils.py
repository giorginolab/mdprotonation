from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import warnings

import MDAnalysis as mda

ResidueKey = tuple[str, int, str]
Coordinate3D = tuple[float, float, float]


@dataclass(frozen=True)
class ParsedPdbAtom:
    residue_key: ResidueKey
    coordinates: Coordinate3D


def parse_pdb_atoms(
    pdb_text: str,
    *,
    filename: str = "structure.pdb",
) -> tuple[ParsedPdbAtom, ...]:
    del filename
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="1 A\\^3 CRYST1 record, this is usually a placeholder.*",
            category=UserWarning,
        )
        universe = mda.Universe(StringIO(pdb_text), format="PDB")
    return tuple(
        ParsedPdbAtom(
            residue_key=_atom_residue_key(atom),
            coordinates=tuple(float(value) for value in atom.position),
        )
        for atom in universe.atoms
    )


def _atom_residue_key(atom: object) -> ResidueKey:
    chain_id = str(getattr(atom, "chainID", "")).strip()
    if not chain_id:
        chain_id = str(getattr(atom, "segid", "")).strip()
    return (
        chain_id or "?",
        int(getattr(atom, "resid")),
        str(getattr(atom, "icode", "")).strip(),
    )
