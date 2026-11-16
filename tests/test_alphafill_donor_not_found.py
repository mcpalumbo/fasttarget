"""An obsolete donor PDB entry (RCSB 404) is a cached final outcome, not an error."""
from ftscripts import alphafill


class Response:
    def __init__(self, status):
        self.status_code = status

    def raise_for_status(self):
        raise AssertionError("404 must be handled before raise_for_status")


def test_rcsb_404_is_not_found_cached_and_not_a_failure(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(alphafill, "_request", lambda session, url, *a: calls.append(url) or Response(404))

    resolver = alphafill.DonorResolver(tmp_path)
    result = resolver.resolve("6g7v", "A")
    assert result["donor_annotation_status"] == "not_found"
    assert len(calls) == 1

    # A later run reads the cached 404 instead of asking RCSB again.
    again = alphafill.DonorResolver(tmp_path).resolve("6G7V", "A")
    assert again["donor_annotation_status"] == "not_found"
    assert len(calls) == 1

    failures = alphafill._execution_failures({}, {"P1": [{"donor_pdb": "6G7V", "donor_protein_chain": "A", **again}]}, {})
    assert failures == []


def test_refresh_cache_asks_again(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(alphafill, "_request", lambda session, url, *a: calls.append(url) or Response(404))
    alphafill.DonorResolver(tmp_path).resolve("6G7V", "A")
    alphafill.DonorResolver(tmp_path, refresh_cache=True).resolve("6G7V", "A")
    assert len(calls) == 2
