from __future__ import annotations

import base64
import hashlib
from dataclasses import dataclass
from math import sqrt

import molviewspec as mvs

from .pdb_utils import ResidueKey, parse_pdb_atoms
from .streamlit_molstar_focus_component import st_molstar_focus_component


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
    del traj_file_content, traj_file_format, traj_file_name
    model_file_name = file_name or f"unknown.{file_format}"
    parser_format = file_format.lower()
    state = _build_scene_state(
        model_file_name=model_file_name,
        parser_format=parser_format,
    )
    assets = {model_file_name: file_content.encode("utf-8")}
    mvsx_bytes = mvs.MVSX(data=state, assets=assets).dumps()
    mvs_data = "base64," + base64.b64encode(mvsx_bytes).decode("ascii")
    model_key = _model_key(
        parser_format=parser_format,
        model_file_name=model_file_name,
        file_content=file_content,
    )
    return st_molstar_focus_component(
        mvs_data=mvs_data,
        model_key=model_key,
        focus_residue=(
            focus_residue.to_component_dict() if focus_residue is not None else None
        ),
        height=_parse_height_px(height),
        key=key,
    )


def _build_scene_state(
    *,
    model_file_name: str,
    parser_format: str,
) -> mvs.State:
    builder = mvs.create_builder()
    structure = (
        builder
        .download(url=model_file_name)
        .parse(format=parser_format)
        .model_structure()
    )
    structure.component(selector="polymer").representation().color(
        custom={"molstar_use_default_coloring": True}
    )
    structure.component(selector="ligand").representation().color(color="blue")
    return builder.get_state()


def _model_key(*, parser_format: str, model_file_name: str, file_content: str) -> str:
    digest = hashlib.sha1(file_content.encode("utf-8"), usedforsecurity=False).hexdigest()
    return f"{parser_format}:{model_file_name}:{digest}"


def _parse_height_px(height: str | int | None) -> int | None:
    if height is None:
        return None
    if isinstance(height, int):
        return height
    value = height.strip()
    if value.endswith("px"):
        value = value[:-2]
    try:
        return int(value)
    except ValueError:
        return None


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
