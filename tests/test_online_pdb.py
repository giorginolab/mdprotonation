from __future__ import annotations

from urllib.error import HTTPError, URLError

import pytest

from pkaScope.online_pdb import OnlinePdbError, fetch_pdb_from_rcsb, normalize_pdb_id


def test_normalize_pdb_id_strips_and_uppercases() -> None:
    assert normalize_pdb_id(" 7bcq ") == "7BCQ"


def test_normalize_pdb_id_rejects_invalid_values() -> None:
    with pytest.raises(ValueError):
        normalize_pdb_id("abc")

    with pytest.raises(ValueError):
        normalize_pdb_id("12345")

    with pytest.raises(ValueError):
        normalize_pdb_id("a!2$")


class _FakeResponse:
    def __init__(self, body: bytes) -> None:
        self._body = body

    def __enter__(self) -> _FakeResponse:
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        del exc_type, exc, tb

    def read(self) -> bytes:
        return self._body


def test_fetch_pdb_from_rcsb_returns_text(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_urlopen(request: object, timeout: float) -> _FakeResponse:
        del timeout
        assert getattr(request, "full_url") == "https://files.rcsb.org/download/7BCQ.pdb"
        return _FakeResponse(b"ATOM      1  N   GLY A   1       0.000   0.000   0.000\n")

    monkeypatch.setattr("pkaScope.online_pdb.urlopen", fake_urlopen)

    pdb_text = fetch_pdb_from_rcsb("7bcq")

    assert "ATOM" in pdb_text


def test_fetch_pdb_from_rcsb_reports_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_urlopen(request: object, timeout: float) -> _FakeResponse:
        del request, timeout
        raise HTTPError(
            url="https://files.rcsb.org/download/XXXX.pdb",
            code=404,
            msg="Not Found",
            hdrs=None,
            fp=None,
        )

    monkeypatch.setattr("pkaScope.online_pdb.urlopen", fake_urlopen)

    with pytest.raises(OnlinePdbError, match="not found"):
        fetch_pdb_from_rcsb("xxxx")


def test_fetch_pdb_from_rcsb_reports_network_errors(monkeypatch: pytest.MonkeyPatch) -> None:
    def fake_urlopen(request: object, timeout: float) -> _FakeResponse:
        del request, timeout
        raise URLError("offline")

    monkeypatch.setattr("pkaScope.online_pdb.urlopen", fake_urlopen)

    with pytest.raises(OnlinePdbError, match="Could not reach RCSB"):
        fetch_pdb_from_rcsb("7BCQ")
