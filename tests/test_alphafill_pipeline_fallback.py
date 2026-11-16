"""Pipeline regression checks for transport cache and one-candidate rescue."""

import json
import gzip
from pathlib import Path

import pandas as pd
import requests

from ftscripts import alphafill


def _blast_row(gene, accession, identity):
    return "\t".join(map(str, [gene, accession, identity, 100, 0, 0,
                               1, 100, 1, 100, 0, 200, 100, 100, 100, 100]))


def test_identity_and_coverage_have_independent_limits(tmp_path):
    blast = tmp_path / "blast.tsv"
    row = _blast_row("G", "VALID", 99).split("\t")
    row[3] = "96"
    row[13] = "96"
    blast.write_text("\t".join(row) + "\n")
    assert alphafill._blast_candidates(blast, {"G"}, 99, 95)["G"][0]["accession"] == "VALID"


def test_pipeline_checks_only_best_homolog_and_reports_absence_separately(tmp_path, monkeypatch):
    files = tmp_path / "Example" / "structures" / "uniprot_files"
    files.mkdir(parents=True)
    (files / "uniprot_species_taxid_1_rest_cdhit100_blast.tsv").write_text(
        "\n".join([_blast_row("G1", "BEST", 100),
                   _blast_row("G1", "NEVER", 99),
                   _blast_row("G2", "MISSING", 98)]) + "\n")
    requested = []

    def fake_download(accession, root, *args):
        requested.append(accession)
        metadata = Path(root) / "cache" / accession / "alphafill.json"
        if accession != "MISSING":
            metadata.parent.mkdir(parents=True, exist_ok=True)
            metadata.write_text(json.dumps({"hits": []}))
        return {"uniprot_id": accession,
                "status": "not_found" if accession == "MISSING" else "no_transplants",
                "metadata_usable": accession != "MISSING",
                "metadata_path": str(metadata) if accession != "MISSING" else "",
                "cif_usable": False, "download_source": "https", "responses": {}}

    monkeypatch.setattr(alphafill, "download_uniprot", fake_download)
    monkeypatch.setattr(alphafill, "_rsync_batch", lambda *args: True)
    alphafill.run_alphafill(str(tmp_path), "Example", {"G1": "DIRECT", "G2": None},
                           {"download_cif": "never", "enrich_donors": False,
                            "generate_pymol": False})
    report = pd.read_csv(tmp_path / "Example" / "ligands" / "alphafill" /
                         "alphafill_selection.tsv", sep="\t", keep_default_na=False)
    by_gene = report.set_index("gene")
    assert by_gene.loc["G1", "selected_uniprot"] == "BEST"
    assert by_gene.loc["G1", "result_status"] == "entry_found_without_transplants"
    assert by_gene.loc["G2", "result_status"] == "alphafill_entry_not_found"
    assert "NEVER" not in requested
    assert by_gene.loc["G1", "candidates_checked"] == 1


def test_valid_archive_json_is_cached_and_avoids_https(tmp_path, monkeypatch):
    root = tmp_path / "alphafill"
    raw = alphafill._archive_path(root, "P12345", "metadata")
    raw.parent.mkdir(parents=True)
    raw.write_text('{"hits": []}')

    def no_http(*args, **kwargs):
        raise AssertionError("cached rsync metadata attempted HTTPS")

    monkeypatch.setattr(requests.Session, "get", no_http)
    first = alphafill.download_uniprot("P12345", root, download_cif="never")
    second = alphafill.download_uniprot("P12345", root, download_cif="never")
    assert first["status"] == second["status"] == "no_transplants"
    assert first["download_source"] == second["download_source"] == "rsync"
    assert second["responses"]["metadata"]["cached"] is True


def test_https_404_is_cached_as_not_found(tmp_path, monkeypatch):
    calls = []

    def missing(self, url, **kwargs):
        calls.append(url)
        response = requests.Response()
        response.status_code = 404
        response.url = url
        return response

    monkeypatch.setattr(requests.Session, "get", missing)
    root = tmp_path / "alphafill"
    first = alphafill.download_uniprot("P12345", root, download_cif="never")
    second = alphafill.download_uniprot("P12345", root, download_cif="never")
    assert first["status"] == second["status"] == "not_found"
    assert len(calls) == 1


def test_rsync_metadata_and_cif_are_reused_on_resume(tmp_path, monkeypatch):
    calls = []

    def fake_rsync(root, accessions, kind, timeout, proxy):
        calls.append((kind, set(accessions)))
        for accession in accessions:
            path = alphafill._archive_path(root, accession, kind)
            path.parent.mkdir(parents=True, exist_ok=True)
            if kind == "metadata":
                path.write_text('{"hits": []}')
            else:
                path.write_bytes(gzip.compress(b"data_test\nloop_\n_atom_site.id\n1\n"))
        return True

    def no_http(*args, **kwargs):
        raise AssertionError("rsync-backed run attempted HTTPS")

    monkeypatch.setattr(alphafill, "_rsync_batch", fake_rsync)
    monkeypatch.setattr(requests.Session, "get", no_http)
    config = {"download_cif": "always", "enrich_donors": False,
              "generate_pymol": False, "similarity_fallback": False}
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert calls == [("metadata", {"P12345"}), ("cif", {"P12345"})]
    calls.clear()
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert calls == []
    root = tmp_path / "Example" / "ligands" / "alphafill"
    assert (root / "models" / "G" / "P12345_alphafill.cif").is_file()
    assert (root / "cache" / "P12345" / "analysis.json").is_file()


def test_successful_ligand_analysis_is_not_repeated(tmp_path, monkeypatch):
    root = tmp_path / "Example" / "ligands" / "alphafill"
    raw = alphafill._archive_path(root, "P12345", "metadata")
    raw.parent.mkdir(parents=True)
    raw.write_text(json.dumps({"hits": [{"pdb_id": "1ABC", "transplants": [{
        "analogue_id": "ATP", "compound_id": "ATP", "asym_id": "A",
        "local_rmsd": 0.5, "clash": {"score": 0.3}}]}]}))
    monkeypatch.setattr(alphafill, "_rsync_batch", lambda *args: True)
    monkeypatch.setattr(requests.Session, "get",
                        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("HTTP requested")))
    original = alphafill.filter_ligand
    calls = []

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(alphafill, "filter_ligand", counted)
    config = {"download_cif": "never", "enrich_donors": False,
              "generate_pymol": False, "similarity_fallback": False}
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert len(calls) == 1
    calls.clear()
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert calls == []


def test_absent_from_successful_rsync_is_not_found_without_https(tmp_path, monkeypatch):
    """rsync --ignore-missing-args skips absent entries: no HTTPS 404 check needed."""
    def no_http(*args, **kwargs):
        raise AssertionError("HTTPS requested for an entry absent from the archive")

    monkeypatch.setattr(alphafill, "_rsync_batch", lambda *args: True)
    monkeypatch.setattr(requests.Session, "get", no_http)
    config = {"download_cif": "always", "enrich_donors": False,
              "generate_pymol": False, "similarity_fallback": False}
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    root = tmp_path / "Example" / "ligands" / "alphafill"
    status = json.loads((root / "cache" / "P12345" / "download_status.json").read_text())
    assert status["status"] == "not_found"
    assert status["download_source"] == "rsync"
    report = pd.read_csv(root / "alphafill_selection.tsv", sep="\t", keep_default_na=False)
    assert report.set_index("gene").loc["G", "result_status"] == "alphafill_entry_not_found"
    # The cached not_found is reused: a second run neither rsyncs nor uses HTTPS.
    monkeypatch.setattr(alphafill, "_rsync_batch",
                        lambda *args: (_ for _ in ()).throw(AssertionError("rsync repeated")))
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)


def test_failed_rsync_still_falls_back_to_https(tmp_path, monkeypatch):
    requested = []

    class Response:
        status_code, url, headers = 404, "", {}

    def fake_get(self, url, **kwargs):
        requested.append(url)
        return Response()

    monkeypatch.setattr(alphafill, "_rsync_batch", lambda *args: False)
    monkeypatch.setattr(requests.Session, "get", fake_get)
    config = {"download_cif": "never", "enrich_donors": False,
              "generate_pymol": False, "similarity_fallback": False}
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert requested and requested[0].endswith("/P12345/json")


def _archive_only(monkeypatch):
    def fake_rsync(root, accessions, kind, timeout, proxy):
        for accession in accessions:
            path = alphafill._archive_path(root, accession, kind)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('{"hits": []}')
        return True
    monkeypatch.setattr(alphafill, "_rsync_batch", fake_rsync)
    monkeypatch.setattr(requests.Session, "get",
                        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("HTTP requested")))


def test_unchanged_completed_run_is_reused_whole(tmp_path, monkeypatch):
    _archive_only(monkeypatch)
    config = {"download_cif": "never", "enrich_donors": False,
              "generate_pymol": False, "similarity_fallback": False}
    first = alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    calls = []
    original = alphafill.download_uniprot
    monkeypatch.setattr(alphafill, "download_uniprot",
                        lambda *args, **kwargs: calls.append(1) or original(*args, **kwargs))
    again = alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert calls == []
    assert again["gene"].tolist() == first["gene"].tolist()

    # Another mapping, or refresh_cache, runs the analysis again.
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345", "H": "Q99999"}, config)
    assert calls
    calls.clear()
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345", "H": "Q99999"},
                            {**config, "refresh_cache": True})
    assert calls


def test_a_run_with_failures_is_not_reused(tmp_path, monkeypatch):
    _archive_only(monkeypatch)
    config = {"download_cif": "never", "enrich_donors": False,
              "generate_pymol": False, "similarity_fallback": False}
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    report_path = tmp_path / "Example" / "ligands" / "alphafill" / "alphafill_results.json"
    report = json.loads(report_path.read_text())
    report["execution"]["status"] = "completed_with_errors"
    report_path.write_text(json.dumps(report))
    calls = []
    original = alphafill.download_uniprot
    monkeypatch.setattr(alphafill, "download_uniprot",
                        lambda *args, **kwargs: calls.append(1) or original(*args, **kwargs))
    alphafill.run_alphafill(str(tmp_path), "Example", {"G": "P12345"}, config)
    assert calls
