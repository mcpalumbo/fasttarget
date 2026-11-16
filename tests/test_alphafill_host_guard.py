import time

import pytest
import requests

from ftscripts import alphafill


class FakeResponse:
    def __init__(self, status=200, headers=None):
        self.status_code, self.headers, self.url = status, headers or {}, ""


class FakeSession:
    def __init__(self, outcomes):
        self.outcomes, self.calls = list(outcomes), 0

    def get(self, url, **kwargs):
        self.calls += 1
        outcome = self.outcomes.pop(0) if self.outcomes else self.default
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    default = FakeResponse()


def refused():
    return requests.ConnectionError("[Errno 111] Connection refused")


def test_breaker_stops_contacting_a_failing_host():
    alphafill._HOSTS.reset(max_failures=3)
    session = FakeSession([refused() for _ in range(10)])
    with pytest.raises(requests.ConnectionError):
        alphafill._request(session, "https://alphafill.eu/v1/aff/P1/json", 1, 5, 0)
    assert session.calls == 3  # stopped at the threshold, not after all 6 attempts
    with pytest.raises(alphafill.HostUnavailable):
        alphafill._request(session, "https://alphafill.eu/v1/aff/P2/json", 1, 5, 0)
    assert session.calls == 3  # skipped without contacting the host
    # Other hosts are not affected.
    assert alphafill._request(FakeSession([]), "https://data.rcsb.org/x", 1, 0, 0).status_code == 200


def test_any_response_resets_the_failure_count():
    alphafill._HOSTS.reset(max_failures=3)
    session = FakeSession([refused(), refused(), FakeResponse(404)])
    assert alphafill._request(session, "https://alphafill.eu/a", 1, 5, 0).status_code == 404
    session = FakeSession([refused(), refused(), FakeResponse()])
    assert alphafill._request(session, "https://alphafill.eu/b", 1, 5, 0).status_code == 200


def test_rate_limit_spaces_requests_to_the_same_host():
    alphafill._HOSTS.reset({"alphafill.eu": 0.2})
    start = time.monotonic()
    for index in range(3):
        alphafill._request(FakeSession([]), f"https://alphafill.eu/{index}", 1, 0, 0)
    assert time.monotonic() - start >= 0.38
    start = time.monotonic()
    alphafill._request(FakeSession([]), "https://data.rcsb.org/x", 1, 0, 0)
    assert time.monotonic() - start < 0.1  # unlimited host


def test_retry_after_is_respected(monkeypatch):
    slept = []
    monkeypatch.setattr(alphafill.time, "sleep", lambda seconds: slept.append(seconds))
    session = FakeSession([FakeResponse(429, {"Retry-After": "7"}), FakeResponse()])
    assert alphafill._request(session, "https://alphafill.eu/x", 1, 2, 0.1).status_code == 200
    assert max(slept) >= 6.9


def test_run_marks_entries_incomplete_when_alphafill_is_unreachable(tmp_path, monkeypatch):
    calls = []

    def refuse(self, url, **kwargs):
        calls.append(url)
        raise refused()

    monkeypatch.setattr(alphafill, "_rsync_batch", lambda *args: False)
    monkeypatch.setattr(requests.Session, "get", refuse)
    mapping = {f"G{index}": f"P{index:05d}" for index in range(40)}
    summary = alphafill.run_alphafill(str(tmp_path), "Example", mapping, {
        "download_cif": "never", "enrich_donors": False, "generate_pymol": False,
        "similarity_fallback": False, "backoff": 0, "https_min_interval": 0,
        "max_connection_failures": 5, "max_workers": 1})
    assert len(calls) == 5  # instead of 40 accessions x 4 attempts
    assert set(summary["alphafill_complete"]) == {False}


@pytest.mark.parametrize("key, value", [("https_min_interval", -1), ("max_connection_failures", 0)])
def test_invalid_guard_options(key, value):
    with pytest.raises(ValueError):
        alphafill.validate_options({key: value})
