from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from pathlib import Path

import streamlit.components.v1 as components
import streamlit_molstar

from .pdb_utils import ResidueKey, parse_pdb_atoms

_COMPONENT_FUNC = components.declare_component(
    "mdprotonation_molstar",
    path=str(Path(streamlit_molstar.__file__).resolve().parent / "frontend" / "build"),
)


@dataclass(frozen=True)
class ResidueFocusTarget:
    chain_id: str
    residue_number: int
    insertion_code: str
    center: tuple[float, float, float]
    radius: float

    @property
    def residue_key(self) -> ResidueKey:
        return (self.chain_id, self.residue_number, self.insertion_code)

    @property
    def token(self) -> str:
        insertion_code = self.insertion_code or "-"
        return f"{self.chain_id}:{self.residue_number}:{insertion_code}"

    def to_component_dict(self) -> dict[str, object]:
        return {
            "token": self.token,
            "chainId": self.chain_id,
            "residueNumber": self.residue_number,
            "insertionCode": self.insertion_code,
            "center": list(self.center),
            "radius": self.radius,
        }


def compute_residue_focus_targets(pdb_text: str) -> dict[ResidueKey, ResidueFocusTarget]:
    residue_atoms: dict[ResidueKey, list[tuple[float, float, float]]] = {}
    for atom in parse_pdb_atoms(pdb_text):
        residue_atoms.setdefault(atom.residue_key, []).append(atom.coordinates)

    focus_targets: dict[ResidueKey, ResidueFocusTarget] = {}
    for residue_key, coordinates in residue_atoms.items():
        center = _centroid(coordinates)
        max_distance = max(_distance(center, point) for point in coordinates)
        focus_targets[residue_key] = ResidueFocusTarget(
            chain_id=residue_key[0],
            residue_number=residue_key[1],
            insertion_code=residue_key[2],
            center=center,
            # Add padding so the camera frames the whole residue comfortably.
            radius=max(2.5, max_distance + 2.0),
        )

    return focus_targets


def st_molstar_focusable_content(
    file_content: str,
    file_format: str,
    traj_file_content: bytes | None = None,
    traj_file_format: str | None = None,
    *,
    file_name: str | None = None,
    traj_file_name: str | None = None,
    focus_residue: ResidueFocusTarget | None = None,
    height: str = "240px",
    key: str | None = None,
):
    params = {
        "scene": "basic",
        "height": height,
        "modelFile": {
            "name": file_name or f"unknown.{file_format}",
            "data": "<placeholder>",
            "format": file_format,
        },
        "modelFile_data": file_content,
        "focusResidue": (
            focus_residue.to_component_dict() if focus_residue is not None else None
        ),
    }
    if traj_file_content is not None and traj_file_format is not None:
        params["trajFile"] = {
            "name": traj_file_name or f"unknown.{traj_file_format}",
            "data": "<placeholder>",
            "format": traj_file_format,
        }
        params["trajFile_data"] = traj_file_content

    return _COMPONENT_FUNC(key=key, default=None, **params)


def _centroid(coordinates: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    count = len(coordinates)
    return tuple(
        sum(point[axis] for point in coordinates) / count for axis in range(3)
    )


def _distance(
    center: tuple[float, float, float],
    point: tuple[float, float, float],
) -> float:
    return sqrt(
        ((point[0] - center[0]) ** 2)
        + ((point[1] - center[1]) ** 2)
        + ((point[2] - center[2]) ** 2)
    )
