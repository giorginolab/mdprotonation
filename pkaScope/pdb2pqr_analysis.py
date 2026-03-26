from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any


@dataclass(frozen=True)
class Pdb2pqrRunOptions:
    ph: float
    force_field: str = "PARSE"
    keep_chain: bool = True
    optimize_hbonds: bool = True
    titration_state_method: str = "propka"


@dataclass(frozen=True)
class Pdb2pqrRunResult:
    optimized_pdb_text: str
    pqr_text: str | None
    atom_count: int
    residue_count: int
    warning_count: int
    diagnostics: tuple[str, ...]


class Pdb2pqrAnalysisError(RuntimeError):
    def __init__(
        self,
        code: str,
        user_message: str,
        *,
        details: str | None = None,
    ) -> None:
        super().__init__(user_message)
        self.code = code
        self.user_message = user_message
        self.details = details


def run_pdb2pqr_analysis(
    pdb_text: str,
    source_name: str,
    options: Pdb2pqrRunOptions,
) -> Pdb2pqrRunResult:
    try:
        run_data = _execute_pdb2pqr_python_api(
            pdb_text=pdb_text,
            source_name=source_name,
            options=options,
        )
    except ModuleNotFoundError as exc:
        raise Pdb2pqrAnalysisError(
            "missing_dependency",
            "pdb2pqr is not installed in this environment.",
            details=str(exc),
        ) from exc
    except Exception as exc:
        raise Pdb2pqrAnalysisError(
            "run_failed",
            "pdb2pqr failed to optimize this structure.",
            details=str(exc),
        ) from exc

    optimized_pdb_text = run_data["optimized_pdb_text"]
    pqr_text = run_data["pqr_text"]
    diagnostics = _split_diagnostics(run_data.get("diagnostics_text", ""))
    atom_count, residue_count = _count_atoms_and_residues(optimized_pdb_text)
    warning_count = sum(
        1
        for line in diagnostics
        if line.startswith("WARNING")
        or "warning" in line.lower()
        or "unable" in line.lower()
        or "missing" in line.lower()
    )
    return Pdb2pqrRunResult(
        optimized_pdb_text=optimized_pdb_text,
        pqr_text=pqr_text,
        atom_count=atom_count,
        residue_count=residue_count,
        warning_count=warning_count,
        diagnostics=diagnostics,
    )


def _execute_pdb2pqr_python_api(
    *,
    pdb_text: str,
    source_name: str,
    options: Pdb2pqrRunOptions,
) -> dict[str, Any]:
    from pdb2pqr.main import run_pdb2pqr

    with TemporaryDirectory(prefix="pkascope-pdb2pqr-") as temp_dir:
        temp_path = Path(temp_dir)
        input_path = temp_path / _safe_name(source_name)
        output_pqr_path = temp_path / f"{input_path.stem}.pqr"
        output_pdb_path = temp_path / f"{input_path.stem}-optimized.pdb"
        input_path.write_text(pdb_text, encoding="utf-8")

        args = [
            "--ff",
            options.force_field,
            "--with-ph",
            f"{options.ph:.2f}",
            "--titration-state-method",
            options.titration_state_method,
            "--pdb-output",
            str(output_pdb_path),
        ]
        if options.keep_chain:
            args.append("--keep-chain")
        if not options.optimize_hbonds:
            args.append("--noopt")
        args.extend([str(input_path), str(output_pqr_path)])

        std_capture = StringIO()
        with redirect_stdout(std_capture), redirect_stderr(std_capture):
            run_pdb2pqr(args)

        if not output_pdb_path.exists():
            raise RuntimeError("pdb2pqr did not produce an optimized PDB output.")
        if not output_pqr_path.exists():
            raise RuntimeError("pdb2pqr did not produce a PQR output.")

        return {
            "optimized_pdb_text": output_pdb_path.read_text(encoding="utf-8"),
            "pqr_text": output_pqr_path.read_text(encoding="utf-8"),
            "diagnostics_text": std_capture.getvalue(),
        }


def _safe_name(source_name: str) -> str:
    base = Path(source_name).name
    if base.lower().endswith(".pdb"):
        return base
    return f"{base}.pdb"


def _split_diagnostics(diagnostics_text: str) -> tuple[str, ...]:
    return tuple(line.strip() for line in diagnostics_text.splitlines() if line.strip())


def _count_atoms_and_residues(pdb_text: str) -> tuple[int, int]:
    atom_count = 0
    residues: set[tuple[str, int, str]] = set()
    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom_count += 1
        chain_id = line[21].strip() or "?"
        try:
            residue_number = int((line[22:26] or "0").strip() or "0")
        except ValueError:
            residue_number = 0
        insertion_code = line[26].strip()
        residues.add((chain_id, residue_number, insertion_code))
    return atom_count, len(residues)
