from __future__ import annotations

import base64
from dataclasses import dataclass
import json
from math import sqrt
from typing import Any

import molviewspec as mvs
import streamlit.components.v1 as components

from .pdb_utils import ResidueKey, parse_pdb_atoms


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
    del traj_file_content, traj_file_format, traj_file_name, key
    model_file_name = file_name or f"unknown.{file_format}"
    parser_format = file_format.lower()
    state = _build_scene_state(
        model_file_name=model_file_name,
        parser_format=parser_format,
        focus_residue=focus_residue,
    )
    assets = {model_file_name: file_content.encode("utf-8")}
    return _render_streamlit_molstar(
        state=state,
        assets=assets,
        height_px=_parse_height_px(height),
        focus_residue=focus_residue,
    )


def _build_scene_state(
    *,
    model_file_name: str,
    parser_format: str,
    focus_residue: ResidueFocusTarget | None,
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
    if focus_residue is not None:
        selector: dict[str, Any] = {
            "auth_seq_id": focus_residue.residue_number,
        }
        if focus_residue.chain_id != "?":
            selector["auth_asym_id"] = focus_residue.chain_id
        if focus_residue.insertion_code:
            selector["pdbx_PDB_ins_code"] = focus_residue.insertion_code
        structure.component(selector=selector).focus(
            radius=max(2.5, focus_residue.radius)
        )
        structure.component(selector=selector).representation(
            type="ball_and_stick",
            size_factor=1.05,
            ignore_hydrogens=True,
        ).color(color="#ff8c00")
    return builder.get_state()


def _render_streamlit_molstar(
    *,
    state: mvs.State,
    assets: dict[str, bytes],
    height_px: int | None,
    focus_residue: ResidueFocusTarget | None,
):
    if focus_residue is None:
        return state.molstar_streamlit(data=assets, height=height_px)

    # Use MolViewSpec for state serialization, then apply camera focus
    # imperatively once the viewer is ready.
    mvsx_bytes = mvs.MVSX(data=state, assets=assets).dumps()
    mvs_data = "base64," + base64.b64encode(mvsx_bytes).decode("ascii")
    focus_payload = {
        "center": [float(v) for v in focus_residue.center],
        "radius": float(max(2.5, focus_residue.radius)),
    }
    html = f"""<!DOCTYPE html>
<html lang="en">
  <head>
    <style>
      #viewer1 {{
        opacity: 0;
      }}
    </style>
    <script src="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"></script>
    <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css" />
  </head>
  <body>
    <div id="viewer1"></div>
    <script>
      const mvsData = {json.dumps(mvs_data)};
      const focus = {json.dumps(focus_payload)};
      molstar.Viewer.create("viewer1", {{
        layoutIsExpanded: false,
        layoutShowControls: false,
        viewportShowToggleFullscreen: true,
        viewportShowExpand: false
      }}).then(async (viewer) => {{
        await viewer.loadMvsData(mvsData, "mvsx");
        window.setTimeout(() => {{
          try {{
            viewer.plugin.managers.camera.focusSphere(
              {{ center: focus.center, radius: focus.radius }},
              {{ durationMs: 0 }}
            );
          }} catch (_err) {{
            // Best effort focus call.
          }}
          const container = document.getElementById("viewer1");
          if (container) container.style.opacity = "1";
        }}, 60);
      }});
    </script>
  </body>
</html>
"""
    return components.html(html, height=height_px)


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
