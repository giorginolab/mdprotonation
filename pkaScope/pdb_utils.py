from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import tempfile
import warnings

import MDAnalysis as mda

ResidueKey = tuple[str, int, str]
Coordinate3D = tuple[float, float, float]


@dataclass(frozen=True)
class ParsedPdbAtom:
    residue_key: ResidueKey
    coordinates: Coordinate3D


def load_pdb_universe(pdb_text: str) -> mda.Universe:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="1 A\\^3 CRYST1 record, this is usually a placeholder.*",
            category=UserWarning,
        )
        return mda.Universe(StringIO(pdb_text), format="PDB")


def write_pdb_atoms(atom_group: object) -> str:
    with tempfile.NamedTemporaryFile(
        mode="w+",
        suffix=".pdb",
        encoding="utf-8",
    ) as handle:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Found no information for attr: 'formalcharges'.*",
                category=UserWarning,
            )
            with mda.Writer(handle.name, format="PDB", multiframe=False) as writer:
                writer.write(atom_group)
        handle.seek(0)
        return handle.read()


def parse_pdb_atoms(
    pdb_text: str,
    *,
    filename: str = "structure.pdb",
) -> tuple[ParsedPdbAtom, ...]:
    del filename
    universe = load_pdb_universe(pdb_text)
    return tuple(
        ParsedPdbAtom(
            residue_key=atom_residue_key(atom),
            coordinates=tuple(float(value) for value in atom.position),
        )
        for atom in universe.atoms
    )


def atom_residue_key(atom: object) -> ResidueKey:
    chain_id = str(getattr(atom, "chainID", "")).strip()
    if not chain_id:
        chain_id = str(getattr(atom, "segid", "")).strip()
    return (
        chain_id or "?",
        int(getattr(atom, "resid")),
        str(getattr(atom, "icode", "")).strip(),
    )
