from __future__ import annotations

from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


class OnlinePdbError(RuntimeError):
    """Raised when a PDB structure cannot be fetched from an online source."""


def normalize_pdb_id(value: str) -> str:
    pdb_id = value.strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError("PDB ID must be exactly 4 letters or digits (for example: 7BCQ).")
    return pdb_id


def fetch_pdb_from_rcsb(pdb_id: str, *, timeout_seconds: float = 10.0) -> str:
    normalized_pdb_id = normalize_pdb_id(pdb_id)
    url = f"https://files.rcsb.org/download/{normalized_pdb_id}.pdb"
    request = Request(url=url, headers={"User-Agent": "mdprotonation/0.1"})

    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            payload = response.read()
    except HTTPError as exc:
        if exc.code == 404:
            raise OnlinePdbError(
                f"PDB ID {normalized_pdb_id} was not found at the RCSB archive."
            ) from exc
        raise OnlinePdbError(
            f"RCSB request failed for {normalized_pdb_id} (HTTP {exc.code})."
        ) from exc
    except URLError as exc:
        reason = str(getattr(exc, "reason", exc))
        raise OnlinePdbError(f"Could not reach RCSB to fetch {normalized_pdb_id}: {reason}") from exc
    except TimeoutError as exc:
        raise OnlinePdbError(f"Timed out while fetching {normalized_pdb_id} from RCSB.") from exc

    try:
        pdb_text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise OnlinePdbError(f"RCSB response for {normalized_pdb_id} was not valid UTF-8.") from exc

    if "ATOM" not in pdb_text and "HETATM" not in pdb_text:
        raise OnlinePdbError(
            f"Downloaded content for {normalized_pdb_id} does not look like a PDB file."
        )
    return pdb_text
