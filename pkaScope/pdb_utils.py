from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
from math import isnan
import tempfile
import warnings

import MDAnalysis as mda

ResidueKey = tuple[str, int, str]
Coordinate3D = tuple[float, float, float]
HISTIDINE_VARIANTS = frozenset({"HID", "HIE", "HIP", "HSD", "HSE", "HSP"})


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


def preprocess_pdb_text(pdb_text: str) -> str:
    universe = load_pdb_universe(pdb_text)
    for residue in universe.residues:
        if str(residue.resname).strip().upper() in HISTIDINE_VARIANTS:
            residue.resname = "HIS"

    selected_atom_indices: dict[
        tuple[str, int, str, str, str],
        tuple[int, float, str],
    ] = {}
    for atom in universe.atoms:
        atom_name = str(getattr(atom, "name", "")).strip().upper()
        if atom_name.startswith("H"):
            continue
        element = str(getattr(atom, "element", "")).strip().upper()
        if element == "H":
            continue

        altloc = str(getattr(atom, "altLoc", "")).strip().upper()
        occupancy = float(getattr(atom, "occupancy", 0.0))
        if isnan(occupancy):
            occupancy = 0.0

        chain_id = str(getattr(atom, "chainID", "")).strip()
        if not chain_id:
            chain_id = str(getattr(atom, "segid", "")).strip()
        key = (
            chain_id or "?",
            int(getattr(atom, "resid")),
            str(getattr(atom, "icode", "")).strip(),
            str(getattr(atom, "resname", "")).strip().upper(),
            atom_name,
        )

        current = selected_atom_indices.get(key)
        if current is None:
            selected_atom_indices[key] = (atom.index, occupancy, altloc)
            continue

        _, current_occupancy, current_altloc = current
        if (
            occupancy > current_occupancy
            or (occupancy == current_occupancy and _altloc_rank(altloc) < _altloc_rank(current_altloc))
        ):
            selected_atom_indices[key] = (atom.index, occupancy, altloc)

    kept_indices = sorted(index for index, _, _ in selected_atom_indices.values())
    return write_pdb_atoms(universe.atoms[kept_indices])


def _altloc_rank(altloc: str) -> tuple[int, str]:
    if altloc == "A":
        return (0, altloc)
    if altloc == "":
        return (1, altloc)
    return (2, altloc)


def atom_residue_key(atom: object) -> ResidueKey:
    chain_id = str(getattr(atom, "chainID", "")).strip()
    if not chain_id:
        chain_id = str(getattr(atom, "segid", "")).strip()
    return (
        chain_id or "?",
        int(getattr(atom, "resid")),
        str(getattr(atom, "icode", "")).strip(),
    )
