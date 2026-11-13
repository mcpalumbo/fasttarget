"""Functions to download, parse, filter and summarize AlphaFill data."""

from __future__ import annotations

import hashlib
import csv
import gzip
import json
import logging
import math
import os
import re
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote

import pandas as pd
import requests


BASE_URL = "https://alphafill.eu/v1/aff"
RSYNC_URL = "rsync://rsync.alphafill.eu/alphafill/"
RCSB_URL = "https://data.rcsb.org/rest/v1/core"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = "FastTarget/AlphaFill"
CONFIDENCE_CATEGORIES = ("high", "intermediate", "low", "unknown")
CONFIDENCE_COUNT_COLUMNS = [f"{category}_confidence_transfers" for category in CONFIDENCE_CATEGORIES]
COMPOUND_SUMMARY_COLUMNS = [
    column for category in CONFIDENCE_CATEGORIES
    for column in (f"alphafill_{category}_confidence_compounds",
                   f"alphafill_{category}_confidence_compound_count")
] + ["alphafill_unique_compounds", "alphafill_unique_accepted_compounds"]
GENE_SUMMARY_COLUMNS = ["gene", "alphafill_status", "alphafill_complete", *COMPOUND_SUMMARY_COLUMNS]
IDENTITY_CATEGORIES = ("below_25", "25-35", "35-45", "45-55", "55-65", "65-75", ">=75", "unknown")
UNIPROT_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,31}$")
LIGAND_COLUMNS = [
    "gene", "uniprot", "selection_source", "mapped_uniprot", "sequence_identity",
    "query_coverage", "subject_coverage", "ligand_instance_id", "compound_code",
    "analogue_id", "donor_compound_id", "asym_id", "chain_id",
    "donor_pdb", "donor_chain", "donor_protein_chain",
    "donor_uniprot", "donor_protein_name", "donor_gene", "donor_organism", "donor_tax_id",
    "donor_annotation_status", "donor_annotation_error", "donor_annotations_json",
    "donor_identity", "global_rmsd", "local_rmsd", "tcs",
    "confidence_category", "identity_category",
    "filter_status", "rejection_reasons", "hit_json", "transplant_json",
    "cif_path", "pymol_script_path",
]
SUMMARY_COLUMNS = [
    "gene", "uniprot", "selection_source", "mapped_uniprot", "sequence_identity",
    "query_coverage", "subject_coverage", "download_source",
    "alphafill_status", "metadata_usable", "total_transfers",
    "accepted_transfers", "rejected_transfers", *CONFIDENCE_COUNT_COLUMNS,
    *COMPOUND_SUMMARY_COLUMNS, "cif_path", "metadata_path", "pymol_script_path",
    "model_output_status", "model_output_error", "error_message",
]


RESULT_PATH_COLUMNS = ("cif_path", "metadata_path", "pymol_script_path")


class AlphaFillExecutionError(RuntimeError):
    """Signal execution failures after available results have been saved."""

    def __init__(self, failures: list[dict[str, Any]], report_path: Path):
        self.failures = failures
        self.report_path = report_path
        examples = "; ".join(f"{item['stage']}: {item['identifier']} ({item['status']})"
                             for item in failures[:5])
        super().__init__(f"AlphaFill encountered {len(failures)} execution failure(s): {examples}. "
                         f"Available results and failure details saved in {report_path}")


def _portable_paths(value: Any, root: Path) -> Any:
    """Serialize local result paths relative to the AlphaFill output directory.

    Runtime paths remain absolute. Only named artifact fields are converted;
    provider URLs, raw metadata and diagnostic messages are left untouched.
    """
    if isinstance(value, dict):
        result = {}
        for key, item in value.items():
            if key in RESULT_PATH_COLUMNS and isinstance(item, str) and item:
                path = Path(item)
                result[key] = path.relative_to(root.resolve()).as_posix() if path.is_absolute() else path.as_posix()
            else:
                result[key] = _portable_paths(item, root)
        return result
    if isinstance(value, list):
        return [_portable_paths(item, root) for item in value]
    return value


def _portable_frame(frame: pd.DataFrame, root: Path) -> pd.DataFrame:
    frame = frame.copy()
    for column in RESULT_PATH_COLUMNS:
        if column in frame:
            frame[column] = frame[column].map(lambda value: _portable_paths({column: value}, root)[column])
    return frame


def _execution_failures(statuses: dict[str, Any], normalized: dict[str, Any],
                        pair_outputs: dict[tuple[str, str], Any]) -> list[dict[str, Any]]:
    """Keep biological absence separate from failed downloads/annotations/exports."""
    failures = []
    for accession, status in sorted(statuses.items()):
        if status.get("status") not in {"available", "no_transplants", "not_found"}:
            failures.append({"stage": "alphafill", "identifier": accession,
                             "status": status.get("status", "unknown_error"),
                             "error": status.get("error_message", "")})
    donors = set()
    for ligands in normalized.values():
        for ligand in ligands:
            state = ligand.get("donor_annotation_status")
            if state == "error" or (state == "partial" and ligand.get("donor_annotation_error")):
                key = (ligand.get("donor_pdb"), ligand.get("donor_protein_chain"))
                if key not in donors:
                    donors.add(key)
                    failures.append({"stage": "donor_annotation", "identifier": "/".join(str(item or "") for item in key),
                                     "status": state, "error": ligand.get("donor_annotation_error", "")})
    for (gene, accession), output in sorted(pair_outputs.items()):
        if output.get("model_output_status") in {"error", "viewer_error"}:
            failures.append({"stage": "model_output", "identifier": f"{gene}/{accession}",
                             "status": output["model_output_status"],
                             "error": output.get("model_output_error", "")})
    return failures


def sanitize_uniprot_id(value: Any) -> str:
    """Return a conservative accession safe for both URLs and paths."""
    accession = str(value or "").strip().upper()
    if not UNIPROT_RE.fullmatch(accession) or ".." in accession:
        raise ValueError(f"Invalid UniProt ID: {value!r}")
    return accession


def _validate_locus_tag(value: Any) -> str:
    """Preserve locus names while rejecting paths outside their model folder."""
    if not isinstance(value, str) or not value.strip() or value in {".", ".."}:
        raise ValueError(f"Invalid locus_tag: {value!r}")
    if any(character in value for character in ("/", "\\", "\0")):
        raise ValueError(f"locus_tag must be a single directory name: {value!r}")
    return value


def summarize_compounds(ligands: list[dict[str, Any]], usable: bool = True) -> dict[str, Any]:
    """Count each placed compound once, at its best known confidence.

    Donor structure, chain, instance and target UniProt do not define a new
    compound. Unknown is used only when no assessed transfer exists. Missing
    results have null counts, distinct from a successful search with no hits.
    """
    priorities = {category: rank for rank, category in enumerate(CONFIDENCE_CATEGORIES)}
    best = {}
    for ligand in ligands:
        code = ligand.get("compound_code")
        if not isinstance(code, str) or not code.strip():
            continue
        code = code.strip().upper()
        category = ligand.get("confidence_category", "unknown")
        if category not in priorities:
            category = "unknown"
        if code not in best or priorities[category] < priorities[best[code]]:
            best[code] = category
    result = {}
    for category in CONFIDENCE_CATEGORIES:
        codes = sorted(code for code, confidence in best.items() if confidence == category)
        result[f"alphafill_{category}_confidence_compounds"] = ";".join(codes)
        result[f"alphafill_{category}_confidence_compound_count"] = len(codes) if usable else None
    result["alphafill_unique_compounds"] = len(best) if usable else None
    result["alphafill_unique_accepted_compounds"] = (
        sum(category in {"high", "intermediate"} for category in best.values()) if usable else None
    )
    return result


def _number(value: Any) -> float | None:
    """
    Convert a value to a finite float when possible.

    :param value: Input value to convert to float.
    :return: Finite float value, or `None` when conversion fails or value is not finite.
    """
    if value is None or isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def normalize_identity(value: Any) -> float | None:
    """Normalize AlphaFill alignment identity to percent (0--100)."""
    identity = _number(value)
    if identity is None or identity < 0 or identity > 100:
        return None
    return identity * 100.0 if identity <= 1.0 else identity


def _canonical(value: Any) -> str:
    """
    Produce a canonical JSON representation for `value`.

    :param value: Any JSON-serializable object.
    :return: Compact JSON string with sorted keys suitable for stable comparisons.
    """
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def parse_alphafill_metadata(uniprot_id: str, payload: Any) -> list[dict[str, Any]]:
    """Translate the AlphaFill JSON schema into one stable model."""
    accession = sanitize_uniprot_id(uniprot_id)
    if not isinstance(payload, dict):
        raise ValueError("AlphaFill metadata must be a JSON object")
    if "hits" not in payload:
        raise ValueError("AlphaFill metadata is missing required field: hits")
    hits = payload["hits"]
    if hits is None:
        return []
    if not isinstance(hits, list):
        raise ValueError("AlphaFill 'hits' must be a list or null")

    candidates = []
    for hit in hits:
        if not isinstance(hit, dict):
            raise ValueError("AlphaFill hit must be an object")
        if "transplants" not in hit:
            raise ValueError("AlphaFill hit is missing required field: transplants")
        transplants = hit["transplants"]
        if transplants is None:
            continue
        if not isinstance(transplants, list):
            raise ValueError("AlphaFill transplants must be a list or null")
        alignment = hit.get("alignment")
        if alignment is not None and not isinstance(alignment, dict):
            raise ValueError("AlphaFill alignment must be an object or null")
        alignment = alignment or {}
        compact_hit = {
            "pdb_id": hit.get("pdb_id"),
            "pdb_asym_id": hit.get("pdb_asym_id"),
            "alignment": alignment,
            "global_rmsd": hit.get("global_rmsd"),
        }
        for transplant in transplants:
            if not isinstance(transplant, dict):
                raise ValueError("AlphaFill transplant must be an object")
            for key in ("analogue_id", "compound_id", "asym_id"):
                value = transplant.get(key)
                if not isinstance(value, str) or not value.strip():
                    raise ValueError(f"AlphaFill transplant requires non-empty string: {key}")
            clash = transplant.get("clash") if isinstance(transplant.get("clash"), dict) else {}
            analogue_id = transplant.get("analogue_id")
            donor_compound_id = transplant.get("compound_id")
            donor_chain = (
                transplant.get("pdb_auth_asym_id")
                or transplant.get("pdb_asym_id")
            )
            # These fields identify the donor residue/instance when supplied by
            # AlphaFill.  Keeping only them avoids coupling the ID to unrelated
            # transplants or future descriptive fields.
            donor_instance = {
                key: transplant.get(key)
                for key in (
                    "id", "pdb_auth_seq_id", "pdb_label_seq_id", "pdb_seq_id",
                    "pdb_auth_ins_code", "pdb_comp_id",
                )
                if transplant.get(key) is not None
            }
            candidates.append({
                "uniprot_id": accession,
                # analogue_id is the compound placed in the AlphaFill model.
                "compound_code": analogue_id or donor_compound_id,
                "analogue_id": analogue_id,
                "donor_compound_id": donor_compound_id,
                "asym_id": transplant.get("asym_id"),
                "chain_id": transplant.get("asym_id"),
                "donor_pdb": hit.get("pdb_id"),
                "donor_chain": donor_chain,
                "donor_protein_chain": hit.get("pdb_asym_id"),
                "donor_instance": donor_instance,
                "donor_identity": normalize_identity(alignment.get("identity")),
                "global_rmsd": _number(hit.get("global_rmsd")),
                "local_rmsd": _number(transplant.get("local_rmsd")),
                "tcs": _number(clash.get("score", transplant.get("tcs"))),
                "hit": compact_hit,
                "transplant": transplant,
            })

    def instance_fingerprint(item: dict[str, Any]) -> str:
        return _canonical({
            "uniprot": accession,
            "donor_pdb": item.get("donor_pdb"),
            "donor_chain": item.get("donor_chain"),
            "analogue_id": item.get("analogue_id"),
            "donor_compound_id": item.get("donor_compound_id"),
            "asym_id": item.get("asym_id"),
            "donor_instance": item.get("donor_instance"),
        })

    # Sorting by the instance fingerprint makes IDs independent of hit and
    # transplant array ordering.
    candidates.sort(key=lambda item: (instance_fingerprint(item), _canonical(item["transplant"])))
    occurrences: dict[str, int] = {}
    for item in candidates:
        fingerprint = instance_fingerprint(item)
        occurrences[fingerprint] = occurrences.get(fingerprint, 0) + 1
        seed = f"{fingerprint}|{occurrences[fingerprint]}"
        item["ligand_instance_id"] = "AF-" + hashlib.sha256(seed.encode()).hexdigest()[:20]
    return candidates


def identity_category(identity: Any) -> str:
    """Bin donor-to-target identity (%), independently of structural confidence.

    Ten-point bins span [25,35) through [65,75), followed by [75,100].
    Values below 25 and missing/invalid identities have separate categories.
    """
    identity = _number(identity)
    if identity is None or not 0 <= identity <= 100:
        return "unknown"
    if identity < 25:
        return "below_25"
    for lower in range(25, 75, 10):
        if identity < lower + 10:
            return f"{lower}-{lower + 10}"
    return ">=75"


def filter_ligand(ligand: dict[str, Any], high_local_rmsd: float = 0.92,
                  high_tcs: float = 0.64, intermediate_local_rmsd: float = 3.10,
                  intermediate_tcs: float = 1.27) -> dict[str, Any]:
    """Categorize structural confidence using strict local RMSD/TCS limits.

    Both metrics must satisfy a tier's limits. Missing or invalid metrics mean
    unknown confidence. Global RMSD and identity never reject a transplant.
    The compatibility field filter_status selects high/intermediate confidence
    for the filtered TSV and accepted_only CIF download policy.
    """
    local_rmsd = _number(ligand.get("local_rmsd"))
    tcs = _number(ligand.get("tcs"))
    reasons = []
    for name, value in (("local_rmsd", local_rmsd), ("tcs", tcs)):
        if value is None:
            reasons.append(f"missing_{name}")
        elif value < 0:
            reasons.append(f"invalid_{name}")
    if reasons:
        confidence = "unknown"
    elif local_rmsd < high_local_rmsd and tcs < high_tcs:
        confidence = "high"
    elif local_rmsd < intermediate_local_rmsd and tcs < intermediate_tcs:
        confidence = "intermediate"
    else:
        confidence = "low"
        if local_rmsd >= intermediate_local_rmsd:
            reasons.append("local_rmsd_at_or_above_intermediate_limit")
        if tcs >= intermediate_tcs:
            reasons.append("tcs_at_or_above_intermediate_limit")
    result = dict(ligand)
    result["confidence_category"] = confidence
    result["identity_category"] = identity_category(ligand.get("donor_identity"))
    result["filter_status"] = "accepted" if confidence in {"high", "intermediate"} else "rejected"
    result["rejection_reasons"] = reasons
    return result


def _atomic_bytes(path: Path, content: bytes) -> None:
    """
    Atomically write `content` bytes to `path`.

    The function writes to a temporary file in the destination directory,
    flushes and fsyncs, then moves the temporary file into place.

    :param path: Destination file path.
    :param content: Bytes to write to the file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _atomic_json(path: Path, value: Any) -> None:
    """
    Atomically write a JSON-serializable `value` to `path` in a human-readable form.

    :param path: Destination file path.
    :param value: JSON-serializable object to write.
    """
    _atomic_bytes(path, (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n").encode())


def _valid_json(path: Path) -> bool:
    """Check that cached JSON has a supported AlphaFill metadata structure."""
    try:
        if not path.is_file() or path.stat().st_size == 0:
            return False
        with path.open(encoding="utf-8") as handle:
            parse_alphafill_metadata(path.parent.name, json.load(handle))
        return True
    except (OSError, ValueError):
        return False


def _valid_cif(path: Path) -> bool:
    """
    Check whether a file has the minimum expected mmCIF structure.

    :param path: Path to the CIF file.
    :return: `True` when it contains a data block and atom-site category.
    """
    try:
        if not path.is_file() or path.stat().st_size == 0:
            return False
        with path.open("rb") as handle:
            prefix = handle.read(1024 * 1024)
        text = prefix.decode("utf-8", errors="replace").lstrip()
        lowered = text[:1024].lower()
        return (
            not lowered.startswith(("<!doctype html", "<html"))
            and "data_" in text
            and "_atom_site." in text
        )
    except OSError:
        return False


def _valid_cif_response(response: requests.Response) -> bool:
    """Validate response headers and content before caching an mmCIF file."""
    content_type = response.headers.get("Content-Type", "").split(";", 1)[0].strip().lower()
    allowed_types = {
        "chemical/x-cif", "chemical/x-mmcif", "text/x-cif",
        "application/x-cif", "text/plain",
        "application/octet-stream",
    }
    if content_type not in allowed_types:
        return False
    text = response.content[:1024 * 1024].decode("utf-8", errors="replace").lstrip()
    lowered = text[:1024].lower()
    return (
        not lowered.startswith(("<!doctype html", "<html"))
        and "data_" in text
        and "_atom_site." in text
    )


def _archive_path(root: Path, accession: str, kind: str) -> Path:
    suffix = "json" if kind == "metadata" else "gz"
    return root / "rsync_raw" / accession[:2] / f"AF-{accession}-F1-filled_v4.cif.{suffix}"


def _source_with(existing: str, source: str) -> str:
    sources = set(existing.split("+")) - {"", "cache"}
    sources.add(source)
    return "+".join(item for item in ("rsync", "https") if item in sources)


def _archive_payload(root: Path, accession: str, kind: str) -> Any | None:
    """Read an already downloaded official archive file, rejecting incomplete files."""
    path = _archive_path(root, accession, kind)
    try:
        if kind == "metadata":
            with path.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            parse_alphafill_metadata(accession, payload)
            return payload
        with gzip.open(path, "rb") as handle:
            data = handle.read()
        if not data.lstrip().startswith(b"data_") or b"_atom_site." not in data:
            return None
        return data
    except (OSError, ValueError, EOFError):
        return None


def _rsync_batch(root: Path, accessions: set[str], kind: str,
                 timeout: float, proxy: str = "") -> bool:
    """Fetch missing archive members in one rsync call; HTTPS handles missing ones."""
    missing = [accession for accession in sorted(accessions)
               if _archive_payload(root, accession, kind) is None]
    if not missing:
        return True
    suffix = "json" if kind == "metadata" else "gz"
    paths = [f"{accession[:2]}/AF-{accession}-F1-filled_v4.cif.{suffix}"
             for accession in missing]
    archive = root / "rsync_raw"
    archive.mkdir(parents=True, exist_ok=True)
    command = ["rsync", "-a", "--ignore-missing-args", "--partial",
               f"--contimeout={max(1, int(timeout))}", f"--timeout={max(1, int(timeout))}",
               "--files-from=-", RSYNC_URL, str(archive) + "/"]
    env = os.environ.copy()
    if proxy:
        if not re.fullmatch(r"[A-Za-z0-9.:-]+", proxy):
            raise ValueError("rsync_proxy must be a host:port SOCKS endpoint")
        env["RSYNC_CONNECT_PROG"] = f"nc -X 5 -x {proxy} %H 873"
    try:
        result = subprocess.run(command, input="\n".join(paths) + "\n", text=True,
                                env=env, capture_output=True, timeout=max(30, timeout * 2 + len(paths)))
    except (OSError, subprocess.TimeoutExpired) as error:
        logging.warning("AlphaFill rsync unavailable: %s", error)
        return False
    if result.returncode:
        logging.warning("AlphaFill rsync failed (%s): %s", result.returncode,
                        result.stderr.strip()[-500:])
        return False
    return True


def _blast_candidates(path: Path, loci: set[str], min_identity: float,
                      min_coverage: float) -> dict[str, list[dict[str, Any]]]:
    """Use full query coverage and aligned length/full subject length."""
    by_gene: dict[str, dict[str, dict[str, Any]]] = {}
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 16 or row[0] not in loci:
                continue
            try:
                accession = sanitize_uniprot_id(row[1])
                identity, aligned, query_coverage, subject_length, bitscore = (
                    float(row[2]), float(row[3]), float(row[13]), float(row[15]), float(row[11]))
                subject_coverage = 100 * aligned / subject_length
            except (ValueError, ZeroDivisionError):
                continue
            if identity < min_identity or min(query_coverage, subject_coverage) < min_coverage:
                continue
            hit = {"accession": accession, "identity": identity,
                   "query_coverage": query_coverage, "subject_coverage": subject_coverage,
                   "bitscore": bitscore}
            old = by_gene.setdefault(row[0], {}).get(accession)
            if old is None or _candidate_score(hit) > _candidate_score(old):
                by_gene[row[0]][accession] = hit
    return {gene: sorted(hits.values(), key=_candidate_score, reverse=True)
            for gene, hits in by_gene.items()}


def _candidate_score(hit: dict[str, Any]) -> tuple[Any, ...]:
    exact = hit["identity"] == hit["query_coverage"] == hit["subject_coverage"] == 100
    return (exact, hit["identity"], min(hit["query_coverage"], hit["subject_coverage"]),
            hit["query_coverage"], hit["subject_coverage"], hit["bitscore"], hit["accession"])


def _request(session: requests.Session, url: str, timeout: float,
             retries: int, backoff: float) -> requests.Response:
    """
    Perform an HTTP GET request with retries and exponential backoff.

    :param session: `requests.Session` instance to use for the request.
    :param url: URL to fetch.
    :param timeout: Timeout in seconds for the request.
    :param retries: Number of retry attempts on transient failures.
    :param backoff: Base backoff interval in seconds (exponential multiplier).
    :return: `requests.Response` on success.
    :raises: Last encountered exception when all retries fail.
    """
    last_error = None
    for attempt in range(retries + 1):
        try:
            response = session.get(url, timeout=timeout, headers={"User-Agent": USER_AGENT})
            if response.status_code != 429 and response.status_code < 500:
                return response
            last_error = requests.HTTPError(f"HTTP {response.status_code}", response=response)
        except (requests.Timeout, requests.ConnectionError) as error:
            last_error = error
        if attempt < retries:
            time.sleep(backoff * (2 ** attempt))
    raise last_error or requests.RequestException("AlphaFill request failed")


def generate_pymol_script(cif_path: str | Path, ligands: list[dict[str, Any]],
                          uniprot_id: str, output_path: str | Path | None = None) -> Path:
    """Write a standalone PyMOL viewer, one independently toggled object per transplant.

    Run the generated file with PyMOL's `run /path/view_alphafill.py` command.
    Ligand asym_id is mmCIF label_asym_id, which PyMOL exposes as segi even
    with its default author-chain numbering. Selecting by chain or compound
    code would incorrectly merge different instances. The original CIF is
    kept intact and loaded as a disabled source object. All classified ligand
    objects start visible; users can toggle each object or confidence group.
    """
    cif_path = Path(cif_path).resolve()
    if not _valid_cif(cif_path):
        raise ValueError(f"Cannot create PyMOL viewer: missing or invalid CIF {cif_path}")
    accession = sanitize_uniprot_id(uniprot_id)
    output_path = Path(output_path).resolve() if output_path else cif_path.with_name("view_alphafill.py")
    relative_cif = os.path.relpath(cif_path, output_path.parent)
    entries = [{key: ligand.get(key) for key in (
        "ligand_instance_id", "compound_code", "asym_id", "confidence_category",
        "donor_pdb", "donor_protein_chain", "donor_protein_name", "donor_organism",
        "donor_identity", "local_rmsd", "tcs",
    )} for ligand in ligands]
    # repr embeds metadata as Python literals, never as executable commands or
    # selections. select_list uses atom indices rather than interpolated IDs.
    source = """# Generated by FastTarget. In PyMOL: run /path/view_alphafill.py
# Toggle individual ligand objects or confidence groups in the right panel.
# Green: high; yellow: intermediate; orange: low; gray: unknown.
def _show_alphafill():
    from pathlib import Path
    import re
    from pymol import cmd

    script_file = globals().get("__script__") or globals().get("__file__")
    if not script_file:
        raise RuntimeError("Open this file with PyMOL's run command")
    cif_path = Path(script_file).resolve().parent / __RELATIVE_CIF__
    if not cif_path.is_file():
        raise FileNotFoundError(cif_path)
    entries = __ENTRIES__
    colors = {"high": "green", "intermediate": "yellow", "low": "orange", "unknown": "gray50"}
    root = cmd.get_unused_name(__ROOT__)
    cmd.group(root)
    source = root + "_source"
    protein = root + "_protein"
    cmd.load(str(cif_path), source, format="cif")
    cmd.hide("everything", source)
    cmd.create(protein, source + " and polymer")
    cmd.show("cartoon", protein)
    cmd.color("lightblue", protein)
    cmd.group(root, protein)
    cmd.group(root, source)
    cmd.disable(source)
    # label_asym_id is always stored as segi, independent of cif_use_auth.
    nonpolymer = {index for _, index in cmd.index(source + " and not polymer")}
    segments = {}
    for atom in cmd.get_model(source).atom:
        if atom.index in nonpolymer:
            segments.setdefault(atom.segi, []).append(atom.index)
    groups = {}
    for category in colors:
        groups[category] = root + "_" + category
        cmd.group(groups[category])
        cmd.group(root, groups[category])
    temporary = root + "_selection"
    for number, entry in enumerate(entries, 1):
        asym_id = str(entry.get("asym_id") or "")
        indices = segments.get(asym_id, [])
        if not indices:
            print("AlphaFill: ligand not found in CIF: " + str(entry.get("ligand_instance_id")) + " / " + asym_id)
            continue
        category = entry.get("confidence_category")
        if category not in colors:
            category = "unknown"
        label = "_".join(str(entry.get(key) or "unknown") for key in
                         ("compound_code", "donor_pdb", "donor_protein_chain", "asym_id"))
        label = re.sub(r"[^A-Za-z0-9_]", "_", label)[:100]
        name = root + "_" + category + "_" + label + "_" + str(number)
        cmd.select_list(temporary, source, indices, mode="index")
        cmd.create(name, temporary)
        cmd.delete(temporary)
        cmd.show("sticks", name)
        cmd.show("spheres", name)
        cmd.set("sphere_scale", 0.25, name)
        cmd.color(colors[category], name)
        cmd.group(groups[category], name)
        print("AlphaFill: " + name + " | " + str(entry))
    cmd.orient(protein)
    print("AlphaFill: toggle each ligand or confidence group in the object panel.")

_show_alphafill()
del _show_alphafill
"""
    source = source.replace("__RELATIVE_CIF__", repr(relative_cif))
    source = source.replace("__ROOT__", repr("AF_" + re.sub(r"[^A-Za-z0-9_]", "_", accession)))
    source = source.replace("__ENTRIES__", repr(entries))
    compile(source, str(output_path), "exec")
    _atomic_bytes(output_path, source.encode("utf-8"))
    return output_path


class DonorResolver:
    """Resolve the donor protein's label_asym_id through RCSB and UniProt.

    Never use the ligand chain or the expression host to identify the protein.
    Cache successful API payloads on disk and all lookup outcomes in memory.
    Multiple UniProt mappings are preserved as separate annotation records.
    """

    def __init__(self, cache_root: str | Path, timeout: float = 30.0,
                 retries: int = 3, backoff: float = 1.0,
                 refresh_cache: bool = False):
        self.cache_root = Path(cache_root)
        self.timeout, self.retries, self.backoff = timeout, retries, backoff
        self.refresh_cache = refresh_cache
        self.session = requests.Session()
        self.payloads: dict[str, Any] = {}
        self.donors: dict[tuple[str, str], dict[str, Any]] = {}

    def close(self) -> None:
        self.session.close()

    def _get(self, url: str, required_key: str, expected_type: type) -> dict[str, Any]:
        if url in self.payloads:
            result = self.payloads[url]
            if isinstance(result, Exception):
                raise result
            return result
        path = self.cache_root / (hashlib.sha256(url.encode()).hexdigest() + ".json")

        def valid(payload: Any) -> bool:
            return isinstance(payload, dict) and isinstance(payload.get(required_key), expected_type)

        try:
            if not self.refresh_cache:
                try:
                    with path.open(encoding="utf-8") as handle:
                        cached = json.load(handle)
                    if valid(cached):
                        self.payloads[url] = cached
                        return cached
                except (OSError, ValueError):
                    pass
            response = _request(self.session, url, self.timeout, self.retries, self.backoff)
            response.raise_for_status()
            payload = response.json()
            if not valid(payload):
                raise ValueError(f"Unsupported donor metadata from {url}: missing/invalid {required_key}")
            _atomic_json(path, payload)
            self.payloads[url] = payload
            return payload
        except (requests.RequestException, OSError, ValueError) as error:
            self.payloads[url] = error
            raise

    def resolve(self, pdb_id: Any, protein_chain: Any) -> dict[str, Any]:
        pdb_id = str(pdb_id or "").strip().upper()
        protein_chain = str(protein_chain or "").strip()
        key = (pdb_id, protein_chain)
        if key in self.donors:
            return self.donors[key]
        result: dict[str, Any] = {
            "donor_uniprot": "", "donor_protein_name": "", "donor_gene": "",
            "donor_organism": "", "donor_tax_id": "",
            "donor_annotation_status": "unresolved", "donor_annotation_error": "",
            "donor_annotations": [],
        }
        self.donors[key] = result
        if not pdb_id or not protein_chain:
            result["donor_annotation_status"] = "missing_identifiers"
            result["donor_annotation_error"] = "Missing donor PDB ID or protein label_asym_id"
            return result
        try:
            instance_url = f"{RCSB_URL}/polymer_entity_instance/{quote(pdb_id, safe='')}/{quote(protein_chain, safe='')}"
            instance = self._get(instance_url, "rcsb_polymer_entity_instance_container_identifiers", dict)
            identifiers = instance["rcsb_polymer_entity_instance_container_identifiers"]
            if identifiers.get("asym_id") != protein_chain:
                raise ValueError("RCSB protein chain does not match donor label_asym_id")
            entity_id = identifiers.get("entity_id")
            if not entity_id:
                raise ValueError("RCSB protein instance has no entity_id")
            entity_url = f"{RCSB_URL}/polymer_entity/{quote(pdb_id, safe='')}/{quote(str(entity_id), safe='')}"
            entity = self._get(entity_url, "rcsb_polymer_entity_container_identifiers", dict)
            entity_identifiers = entity["rcsb_polymer_entity_container_identifiers"]
            # RCSB may omit uniprot_ids when there is no mapping. Keep all
            # mapped accessions, including chimeric/multi-protein entities.
            accessions = entity_identifiers.get("uniprot_ids") or []
            if not isinstance(accessions, list):
                raise ValueError("RCSB uniprot_ids must be a list")
            if not accessions:
                accessions = [item["database_accession"]
                              for item in entity_identifiers.get("reference_sequence_identifiers", [])
                              if item.get("database_name") == "UniProt" and item.get("database_accession")]
            accessions = sorted({sanitize_uniprot_id(value) for value in accessions})
            result["donor_uniprot"] = ";".join(accessions)
            if not accessions:
                result["donor_annotation_status"] = "unmapped"
                result["donor_annotation_error"] = "No UniProt mapping for the donor protein entity"
                return result
            records, errors = [], []
            for accession in accessions:
                url = f"{UNIPROT_URL}/{quote(accession, safe='')}.json"
                record = {"uniprot": accession, "protein_name": "", "gene": "",
                          "organism": "", "tax_id": "", "status": "error", "error": ""}
                try:
                    data = self._get(url, "primaryAccession", str)
                    description = data.get("proteinDescription") or {}
                    name = description.get("recommendedName") or next(iter(description.get("submissionNames") or []), {})
                    record["protein_name"] = (name.get("fullName") or {}).get("value", "")
                    genes = []
                    for gene in data.get("genes") or []:
                        if (gene.get("geneName") or {}).get("value"):
                            genes.append(gene["geneName"]["value"])
                        else:
                            genes.extend(item["value"] for field in ("orderedLocusNames", "orfNames")
                                         for item in gene.get(field) or [] if item.get("value"))
                    record["gene"] = ";".join(sorted(set(genes)))
                    organism = data.get("organism") or {}
                    record["organism"] = organism.get("scientificName", "")
                    record["tax_id"] = str(organism.get("taxonId") or "")
                    record["resolved_uniprot"] = data["primaryAccession"]
                    record["status"] = "available" if record["protein_name"] and record["organism"] else "partial"
                except (requests.RequestException, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
                    record["error"] = str(error)
                    errors.append(f"{accession}: {error}")
                records.append(record)
            result["donor_annotations"] = records
            for output, field in (("donor_protein_name", "protein_name"), ("donor_gene", "gene"),
                                  ("donor_organism", "organism"), ("donor_tax_id", "tax_id")):
                result[output] = ";".join(sorted({record[field] for record in records if record[field]}))
            states = {record["status"] for record in records}
            result["donor_annotation_status"] = ("available" if states == {"available"}
                                                  else "error" if states == {"error"} else "partial")
            result["donor_annotation_error"] = "; ".join(errors)
        except (requests.RequestException, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
            result["donor_annotation_status"] = "error"
            result["donor_annotation_error"] = str(error)
            logging.warning("Donor lookup failed for %s/%s: %s", pdb_id, protein_chain, error)
        return result


def download_uniprot(uniprot_id: str, output_root: str | Path,
                     timeout: float = 30.0, retries: int = 3,
                     backoff: float = 1.0, download_cif: str = "always",
                     thresholds: dict[str, float] | None = None,
                     refresh_cache: bool = False,
                     locus_tag: str | None = None) -> dict[str, Any]:
    """Fetch one accession; reuse valid cache unless refresh_cache is True.

    output_root contains models/<locus_tag>/<accession>_alphafill.cif as a primary result
    and cache/<accession>/ metadata/status files.
    A failed refresh preserves previous files on disk, but marks them unusable
    for this run. CIF filtering selects downloads, never edits their contents.
    """
    accession = sanitize_uniprot_id(uniprot_id)
    if download_cif not in {"always", "accepted_only", "never"}:
        raise ValueError(f"Invalid download_cif policy: {download_cif!r}")
    if not isinstance(refresh_cache, bool):
        raise ValueError("refresh_cache must be a boolean")
    output_root = Path(output_root).resolve()
    directory = output_root / "cache" / accession
    model_directory = Path(output_root) / "models" / _validate_locus_tag(locus_tag if locus_tag is not None else accession)
    cif_path = model_directory / f"{accession}_alphafill.cif"
    metadata_path = directory / "alphafill.json"
    status_path = directory / "download_status.json"
    directory.mkdir(parents=True, exist_ok=True)
    urls = {
        "cif": f"{BASE_URL}/{quote(accession, safe='')}",
        "metadata": f"{BASE_URL}/{quote(accession, safe='')}/json",
    }
    status: dict[str, Any] = {
        "uniprot_id": accession, "status": "error_transient", "urls": urls,
        "cif_path": "", "metadata_path": "",
        "checked_at": datetime.now(timezone.utc).isoformat(),
        "responses": {}, "error_message": "", "refresh_cache": refresh_cache,
        "metadata_usable": False, "cif_usable": False, "path_base": "../..",
        "download_source": "",
    }
    session = requests.Session()

    def fetch(kind: str) -> requests.Response:
        response = _request(session, urls[kind], timeout, retries, backoff)
        status["responses"][kind] = {
            "url": response.url, "status_code": response.status_code,
            "content_type": response.headers.get("Content-Type"),
            "etag": response.headers.get("ETag"),
            "last_modified": response.headers.get("Last-Modified"), "cached": False,
        }
        return response

    try:
        previous = {}
        if not refresh_cache and status_path.is_file():
            try:
                previous = json.loads(status_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                pass
        if not refresh_cache and _valid_json(metadata_path):
            with metadata_path.open(encoding="utf-8") as handle:
                payload = json.load(handle)
            status["responses"]["metadata"] = {"url": urls["metadata"], "cached": True}
            status["download_source"] = previous.get("download_source", "cache")
        elif (not refresh_cache and previous.get("status") == "not_found"
              and _archive_payload(output_root, accession, "metadata") is None):
            status["status"] = "not_found"
            status["responses"]["metadata"] = {"cached": True, "status_code": 404}
            status["error_message"] = "AlphaFill metadata unavailable (cached HTTP 404)"
            return status
        else:
            payload = None if refresh_cache else _archive_payload(output_root, accession, "metadata")
            if payload is not None:
                status["responses"]["metadata"] = {"source": RSYNC_URL, "cached": True}
                status["download_source"] = "rsync"
            else:
                response = fetch("metadata")
                status["download_source"] = "https"
                if response.status_code == 404:
                    status["status"] = "not_found"
                    status["error_message"] = "AlphaFill metadata unavailable (HTTP 404)"
                    return status
                response.raise_for_status()
                try:
                    payload = response.json()
                except ValueError as error:
                    status["status"] = "invalid_json"
                    status["error_message"] = str(error)
                    return status
                try:
                    parse_alphafill_metadata(accession, payload)
                except ValueError as error:
                    status["status"] = "unsupported_schema"
                    status["error_message"] = str(error)
                    return status
            _atomic_json(metadata_path, payload)

        ligands = parse_alphafill_metadata(accession, payload)
        status["metadata_usable"] = True
        status["metadata_path"] = str(metadata_path)
        cif_required = download_cif == "always" or (
            download_cif == "accepted_only"
            and any(filter_ligand(item, **(thresholds or {}))["filter_status"] == "accepted"
                    for item in ligands)
        )
        if not cif_required:
            status["responses"]["cif"] = {
                "url": urls["cif"], "cached": False,
                "skipped": True, "policy": download_cif,
            }
        elif not refresh_cache and _valid_cif(cif_path):
            status["responses"]["cif"] = {"url": urls["cif"], "cached": True}
            status["cif_usable"] = True
            status["cif_path"] = str(cif_path)
        else:
            raw_cif = None if refresh_cache else _archive_payload(output_root, accession, "cif")
            if raw_cif is not None:
                status["responses"]["cif"] = {"source": RSYNC_URL, "cached": True}
                status["download_source"] = _source_with(status["download_source"], "rsync")
            else:
                response = fetch("cif")
                if response.status_code == 404:
                    status["status"] = "cif_missing"
                    status["error_message"] = "AlphaFill metadata exists but CIF is unavailable (HTTP 404)"
                    return status
                response.raise_for_status()
                if not _valid_cif_response(response):
                    status["status"] = "invalid_cif"
                    status["error_message"] = "AlphaFill CIF response failed Content-Type/mmCIF validation"
                    return status
                raw_cif = response.content
                status["download_source"] = _source_with(status["download_source"], "https")
            _atomic_bytes(cif_path, raw_cif)
            status["cif_usable"] = True
            status["cif_path"] = str(cif_path)
        status["status"] = "available" if ligands else "no_transplants"
    except requests.RequestException as error:
        response = error.response
        if response is not None and 400 <= response.status_code < 500 and response.status_code != 429:
            status["status"] = "http_error"
        status["error_message"] = str(error)
        logging.warning("AlphaFill failed for %s: %s", accession, error)
    except (OSError, ValueError) as error:
        status["status"] = "invalid_metadata" if isinstance(error, ValueError) else "io_error"
        status["error_message"] = str(error)
        logging.warning("AlphaFill failed for %s: %s", accession, error)
    finally:
        session.close()
        _atomic_json(status_path, _portable_paths(status, output_root))
    return status


def _mapping_pairs(mapping: dict[str, Any]) -> list[tuple[str, str]]:
    """
    Convert a gene->accession mapping into a sorted list of (gene, accession) pairs.

    :param mapping: Mapping where keys are gene names and values are accession or list of accessions.
    :return: Sorted unique list of `(gene, accession)` pairs with sanitized UniProt accessions.
    """
    pairs = []
    for gene in sorted(mapping):
        values = mapping[gene]
        values = values if isinstance(values, list) else [values]
        for value in values:
            if value:
                try:
                    pairs.append((str(gene), sanitize_uniprot_id(value)))
                except ValueError as error:
                    logging.warning("Skipping AlphaFill mapping %s: %s", gene, error)
    return sorted(set(pairs))


def _tsv_row(gene: str, accession: str, ligand: dict[str, Any]) -> dict[str, Any]:
    """
    Format a ligand dictionary into a row suitable for TSV output.

    :param gene: Gene name associated with the ligand.
    :param accession: UniProt accession for the gene.
    :param ligand: Ligand dictionary produced by `parse_alphafill_metadata`/`filter_ligand`.
    :return: Dictionary with keys matching `LIGAND_COLUMNS` for CSV/TSV writing.
    """
    return {
        "gene": gene, "uniprot": accession,
        **{key: ligand.get(key) for key in LIGAND_COLUMNS if key not in {"gene", "uniprot", "rejection_reasons", "hit_json", "transplant_json", "donor_annotations_json"}},
        "rejection_reasons": ";".join(ligand.get("rejection_reasons", [])),
        "donor_annotations_json": _canonical(ligand.get("donor_annotations", [])),
        "hit_json": _canonical(ligand.get("hit")),
        "transplant_json": _canonical(ligand.get("transplant")),
    }


def _validated_config(config: dict[str, Any]) -> tuple[dict[str, float], dict[str, Any]]:
    """Convert and validate AlphaFill configuration in one place."""
    obsolete = {"min_identity", "max_global_rmsd", "max_local_rmsd"}.intersection(config)
    if obsolete:
        raise ValueError("AlphaFill now uses confidence categories; remove obsolete filters: "
                         + ", ".join(sorted(obsolete)))
    if "retry_backoff" in config:
        raise ValueError("Use AlphaFill config key 'backoff', not 'retry_backoff'")
    try:
        thresholds = {
            "high_local_rmsd": float(config.get("high_local_rmsd", 0.92)),
            "high_tcs": float(config.get("high_tcs", 0.64)),
            "intermediate_local_rmsd": float(config.get("intermediate_local_rmsd", 3.10)),
            "intermediate_tcs": float(config.get("intermediate_tcs", 1.27)),
        }
        options = {
            "timeout": float(config.get("timeout", 30.0)),
            "retries": int(config.get("retries", 3)),
            "backoff": float(config.get("backoff", 1.0)),
            "max_workers": int(config.get("max_workers", 4)),
            "download_cif": str(config.get("download_cif", "accepted_only")),
            "refresh_cache": config.get("refresh_cache", False),
            "enrich_donors": config.get("enrich_donors", True),
            "generate_pymol": config.get("generate_pymol", True),
            "strict_errors": config.get("strict_errors", False),
            "rsync_proxy": str(config.get("rsync_proxy", "")),
            "fallback_identity": float(config.get("fallback_identity", 95.0)),
            "fallback_coverage": float(config.get("fallback_coverage", 95.0)),
            "similarity_fallback": config.get("similarity_fallback", True),
        }
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid AlphaFill configuration value: {error}") from error
    errors = []
    if not isinstance(options["refresh_cache"], bool):
        errors.append("refresh_cache must be a boolean")
    if not isinstance(options["enrich_donors"], bool):
        errors.append("enrich_donors must be a boolean")
    if not isinstance(options["generate_pymol"], bool):
        errors.append("generate_pymol must be a boolean")
    if not isinstance(options["strict_errors"], bool):
        errors.append("strict_errors must be a boolean")
    if not isinstance(options["similarity_fallback"], bool):
        errors.append("similarity_fallback must be a boolean")
    if options["rsync_proxy"] and not re.fullmatch(r"[A-Za-z0-9.:-]+", options["rsync_proxy"]):
        errors.append("rsync_proxy must be a host:port SOCKS endpoint")
    for key in ("fallback_identity", "fallback_coverage"):
        if not math.isfinite(options[key]) or not 0 < options[key] <= 100:
            errors.append(f"{key} must be between 0 and 100")
    for key, value in thresholds.items():
        if not math.isfinite(value) or value <= 0:
            errors.append(f"{key} must be a positive finite number")
    for metric in ("local_rmsd", "tcs"):
        if thresholds[f"high_{metric}"] >= thresholds[f"intermediate_{metric}"]:
            errors.append(f"high_{metric} must be less than intermediate_{metric}")
    if not math.isfinite(options["timeout"]) or options["timeout"] <= 0:
        errors.append("timeout must be a positive finite number")
    if options["retries"] < 0:
        errors.append("retries must be >= 0")
    if not math.isfinite(options["backoff"]) or options["backoff"] < 0:
        errors.append("backoff must be a non-negative finite number")
    if not 1 <= options["max_workers"] <= 64:
        errors.append("max_workers must be between 1 and 64")
    if options["download_cif"] not in {"always", "accepted_only", "never"}:
        errors.append("download_cif must be one of: always, accepted_only, never")
    if errors:
        raise ValueError("Invalid AlphaFill configuration: " + "; ".join(errors))
    return thresholds, options


def validate_options(config: dict[str, Any]) -> dict[str, Any]:
    """Validate and normalize options for callers such as the YAML loader."""
    thresholds, options = _validated_config(config)
    return {**thresholds, **options}


def run_alphafill(output_path: str, organism_name: str, mapping: dict[str, Any],
                  config: dict[str, Any]) -> pd.DataFrame:
    """Write organism outputs; config refresh_cache=True bypasses cached downloads.

    Genes without a usable accession are retained with status unmapped.
    enrich_donors (default True) resolves donor protein chains via RCSB/UniProt.
    Raw donor API responses are cached in donor_cache and refreshed with refresh_cache.
    CIF models are primary outputs under ligands/alphafill/models/<locus_tag>/.
    Return one row per locus with unique compound lists/counts; accession-level
    transfer counts and model paths are written to alphafill_accession_summary.tsv.
    generate_pymol (default True) writes a viewer beside each usable CIF; its path
    is included in the accession summary. download_cif=always also covers low-only entries.
    Serialized paths are relative to ligands/alphafill/ (not the TSV/JSON file
    location). strict_errors=True raises AlphaFillExecutionError after writing
    outputs if downloads, donor annotation or exports failed; absent AlphaFill
    entries and missing donor mappings remain valid outcomes. Default: False.
    """
    thresholds, options = _validated_config(config)
    for gene in mapping:
        _validate_locus_tag(gene)
    root = Path(output_path).resolve() / organism_name / "ligands" / "alphafill"
    cache_root = root / "cache"
    root.mkdir(parents=True, exist_ok=True)
    pairs = _mapping_pairs(mapping)
    direct_pairs = list(pairs)
    accessions = sorted({accession for _, accession in pairs})
    primary_loci = {}
    for gene, accession in pairs:
        primary_loci.setdefault(accession, gene)
    statuses = {}
    rsync_available = True

    def download_batch(wanted: set[str], policy: str) -> None:
        nonlocal rsync_available
        if not wanted:
            return
        metadata_needed = set()
        for accession in wanted:
            cached = cache_root / accession
            if not options["refresh_cache"] and _valid_json(cached / "alphafill.json"):
                continue
            try:
                previous = json.loads((cached / "download_status.json").read_text(encoding="utf-8"))
            except (OSError, ValueError):
                previous = {}
            if not options["refresh_cache"] and previous.get("status") == "not_found":
                continue
            metadata_needed.add(accession)
        if rsync_available and metadata_needed:
            rsync_available = _rsync_batch(root, metadata_needed, "metadata",
                                           options["timeout"], options["rsync_proxy"])
        def execute(accessions_for_call: set[str], policy_for_call: str, refresh: bool) -> None:
            if not accessions_for_call:
                return
            with ThreadPoolExecutor(max_workers=min(options["max_workers"], len(accessions_for_call))) as executor:
                futures = {executor.submit(download_uniprot, accession, root,
                                           options["timeout"], options["retries"],
                                           options["backoff"], policy_for_call,
                                           thresholds, refresh, primary_loci[accession]): accession
                           for accession in accessions_for_call}
                for future in as_completed(futures):
                    accession = futures[future]
                    try:
                        statuses[accession] = future.result()
                    except Exception as error:
                        logging.exception("Unexpected AlphaFill error for %s", accession)
                        statuses[accession] = {"uniprot_id": accession, "status": "error_transient",
                                               "error_message": str(error), "responses": {}}

        execute(wanted, "never", options["refresh_cache"])
        if policy == "never":
            return
        usable = {accession for accession in wanted if statuses[accession].get("metadata_usable")}
        if rsync_available and policy != "never":
            cif_needed = set()
            for accession in usable:
                if _valid_cif(root / "models" / primary_loci[accession] / f"{accession}_alphafill.cif") and not options["refresh_cache"]:
                    continue
                metadata = cache_root / accession / "alphafill.json"
                payload = None
                try:
                    payload = json.loads(metadata.read_text(encoding="utf-8")) if _valid_json(metadata) else None
                    ligands = parse_alphafill_metadata(accession, payload) if payload is not None else []
                except (OSError, ValueError):
                    ligands = []
                if payload is not None and (policy == "always" or any(
                        filter_ligand(item, **thresholds)["filter_status"] == "accepted" for item in ligands)):
                    cif_needed.add(accession)
            if cif_needed:
                rsync_available = _rsync_batch(root, cif_needed, "cif",
                                               options["timeout"], options["rsync_proxy"])
        execute(usable, policy, False)

    download_batch(set(accessions), options["download_cif"])

    selection: dict[str, dict[str, Any]] = {}
    for gene, accession in direct_pairs:
        status = statuses[accession]
        if not status.get("metadata_usable"):
            continue
        try:
            with (cache_root / accession / "alphafill.json").open(encoding="utf-8") as handle:
                count = len(parse_alphafill_metadata(accession, json.load(handle)))
        except (OSError, ValueError):
            count = 0
        if count and gene not in selection:
            selection[gene] = {"accession": accession, "source": "direct", "identity": "",
                               "query_coverage": "", "subject_coverage": "",
                               "transplants": count}

    unresolved = set(mapping) - set(selection)
    blast_files = sorted((Path(output_path).resolve() / organism_name / "structures" / "uniprot_files").glob(
        "uniprot_species_taxid_*_rest_cdhit100_blast.tsv"))
    blast_path = blast_files[0] if blast_files else None
    candidates = (_blast_candidates(blast_path, unresolved, options["fallback_identity"],
                                    options["fallback_coverage"])
                  if options["similarity_fallback"] and blast_path else {})
    direct_by_gene = {gene: {accession for g, accession in direct_pairs if g == gene}
                      for gene in mapping}
    for gene in candidates:
        candidates[gene] = [hit for hit in candidates[gene]
                            if hit["accession"] not in direct_by_gene[gene]]
    checked = {gene: 0 for gene in mapping}
    pending = {gene: hits[0] for gene, hits in candidates.items()
               if gene not in selection and hits}
    for gene, hit in pending.items():
        primary_loci.setdefault(hit["accession"], gene)
    download_batch({hit["accession"] for hit in pending.values()} - set(statuses), "never")
    for gene, hit in pending.items():
        accession = hit["accession"]
        checked[gene] = 1
        count = 0
        if statuses[accession].get("metadata_usable"):
            try:
                with (cache_root / accession / "alphafill.json").open(encoding="utf-8") as handle:
                    count = len(parse_alphafill_metadata(accession, json.load(handle)))
            except (OSError, ValueError):
                pass
        exact = hit["identity"] == hit["query_coverage"] == hit["subject_coverage"] == 100
        selection[gene] = {**hit, "source": "exact" if exact else "similar",
                           "transplants": count}
    fallback_accessions = {hit["accession"] for hit in selection.values() if hit["source"] != "direct"}
    if options["download_cif"] != "never":
        download_batch(fallback_accessions, options["download_cif"])
    pairs = sorted(set(pairs) | {(gene, hit["accession"]) for gene, hit in selection.items()
                                 if hit["source"] != "direct"})
    accessions = sorted({accession for _, accession in pairs})

    selection_rows = []
    for gene in sorted(mapping):
        hit = selection.get(gene, {})
        matched = hit.get("accession", "")
        states = [statuses[accession].get("status", "error_transient")
                  for accession in direct_by_gene[gene]]
        if matched and hit.get("transplants"):
            result_status = "transplants_found"
        elif matched and statuses[matched].get("status") == "not_found":
            result_status = "alphafill_entry_not_found"
        elif matched and statuses[matched].get("status") not in {"available", "no_transplants", "not_found"}:
            result_status = "incomplete_download"
        elif matched:
            result_status = "entry_found_without_transplants"
        elif any(state not in {"available", "no_transplants", "not_found"} for state in states):
            result_status = "incomplete_download"
        elif states and all(state == "not_found" for state in states):
            result_status = "alphafill_entry_not_found"
        elif not direct_by_gene[gene] and not candidates.get(gene):
            result_status = "unmapped_or_no_candidates"
        else:
            result_status = "entry_found_without_transplants"
        selection_rows.append({"gene": gene, "mapped_uniprot": ";".join(sorted(direct_by_gene[gene])),
                               "selected_uniprot": matched, "selection_source": hit.get("source", "none"),
                               "mapped_alphafill_status": ";".join(sorted(set(states))),
                               "selected_alphafill_status": statuses.get(matched, {}).get("status", ""),
                               "identity": hit.get("identity", ""),
                               "query_coverage": hit.get("query_coverage", ""),
                               "subject_coverage": hit.get("subject_coverage", ""),
                               "transplants": hit.get("transplants", ""),
                               "candidates_checked": checked[gene], "result_status": result_status,
                               "download_source": statuses.get(matched, {}).get("download_source", "")})
    _atomic_bytes(root / "alphafill_selection.tsv", pd.DataFrame(selection_rows).to_csv(
        sep="\t", index=False).encode())

    normalized: dict[str, list[dict[str, Any]]] = {}
    analysis_keys: dict[str, str] = {}
    restored_analysis: set[str] = set()
    for accession in accessions:
        metadata_path = cache_root / accession / "alphafill.json"
        status = statuses[accession]
        if not status.get("metadata_usable", False):
            normalized[accession] = []
            continue
        try:
            cache_key = hashlib.sha256(metadata_path.read_bytes() + json.dumps(
                {"thresholds": thresholds, "enrich_donors": options["enrich_donors"]},
                sort_keys=True).encode()).hexdigest()
            analysis_keys[accession] = cache_key
            analysis_path = cache_root / accession / "analysis.json"
            if not options["refresh_cache"] and analysis_path.is_file():
                try:
                    cached_analysis = json.loads(analysis_path.read_text(encoding="utf-8"))
                    if cached_analysis.get("key") == cache_key and isinstance(cached_analysis.get("ligands"), list):
                        normalized[accession] = cached_analysis["ligands"]
                        restored_analysis.add(accession)
                        continue
                except (OSError, ValueError):
                    pass
            with metadata_path.open(encoding="utf-8") as handle:
                parsed = parse_alphafill_metadata(accession, json.load(handle))
            normalized[accession] = [filter_ligand(item, **thresholds) for item in parsed]
        except (OSError, ValueError, json.JSONDecodeError) as error:
            normalized[accession] = []
            status = statuses[accession]
            if isinstance(error, json.JSONDecodeError):
                status["status"] = "parse_error"
            elif isinstance(error, OSError):
                status["status"] = (
                    status.get("status")
                    if status.get("status") not in {None, "available"}
                    else "invalid_metadata"
                )
            else:
                status["status"] = "unsupported_schema"
            status["error_message"] = str(error)
            status["metadata_usable"] = False
            status["metadata_path"] = ""

    # Resolve each unique donor protein independently of ligand quality.
    # Annotation failures must not change AlphaFill classifications or counts.
    resolver = DonorResolver(root / "donor_cache", options["timeout"],
                             options["retries"], options["backoff"], options["refresh_cache"])
    try:
        for accession, ligands in normalized.items():
            if accession in restored_analysis:
                continue
            for ligand in ligands:
                if options["enrich_donors"]:
                    ligand.update(resolver.resolve(ligand.get("donor_pdb"), ligand.get("donor_protein_chain")))
                else:
                    ligand["donor_annotation_status"] = "not_requested"
            if accession in analysis_keys and all(
                    ligand.get("donor_annotation_status") != "error" and
                    not (ligand.get("donor_annotation_status") == "partial" and ligand.get("donor_annotation_error"))
                    for ligand in ligands):
                _atomic_json(cache_root / accession / "analysis.json",
                             {"key": analysis_keys[accession], "ligands": ligands})
    finally:
        resolver.close()

    # A download is shared by accession, but outputs and viewers belong to loci.
    pair_outputs = {}
    for gene, accession in pairs:
        status = statuses[accession]
        output = {"cif_path": "", "pymol_script_path": "",
                  "model_output_status": "not_available", "model_output_error": ""}
        if status.get("cif_usable") and status.get("metadata_usable"):
            destination = root / "models" / gene / f"{accession}_alphafill.cif"
            try:
                source = Path(status["cif_path"])
                if source.resolve() != destination.resolve():
                    _atomic_bytes(destination, source.read_bytes())
                output["cif_path"] = str(destination)
                output["model_output_status"] = "available"
            except OSError as error:
                output["model_output_status"] = "error"
                output["model_output_error"] = str(error)
            if output["cif_path"] and options["generate_pymol"]:
                try:
                    viewer = generate_pymol_script(
                        destination, normalized[accession], accession,
                        destination.with_name(f"{accession}_view_alphafill.py"),
                    )
                    output["pymol_script_path"] = str(viewer)
                except (OSError, ValueError) as error:
                    output["model_output_status"] = "viewer_error"
                    output["model_output_error"] = str(error)
        pair_outputs[(gene, accession)] = output
        status.setdefault("model_outputs", {})[gene] = output

    for accession, status in statuses.items():
        status["path_base"] = "../.."
        _atomic_json(cache_root / accession / "download_status.json", _portable_paths(status, root))

    all_rows, summary_rows, genes_json = [], [], {}
    mapped_genes = {gene for gene, _ in pairs}
    for gene in sorted(set(mapping) - mapped_genes):
        summary_rows.append({
            "gene": gene, "uniprot": "", "selection_source": "none", "mapped_uniprot": "",
            "sequence_identity": "", "query_coverage": "", "subject_coverage": "",
            "download_source": "", "alphafill_status": "unmapped", "metadata_usable": False,
            "total_transfers": 0, "accepted_transfers": 0, "rejected_transfers": 0,
            **{column: 0 for column in CONFIDENCE_COUNT_COLUMNS},
            **summarize_compounds([], usable=False),
            "cif_path": "", "metadata_path": "", "pymol_script_path": "",
            "model_output_status": "not_available", "model_output_error": "",
            "error_message": "No valid UniProt accession in gene mapping",
        })
        genes_json[gene] = {"": {
            "download": {"status": "unmapped", "error_message": "No valid UniProt accession in gene mapping"},
            "thresholds": thresholds, "ligands": [],
        }}
    for gene, accession in pairs:
        ligands = normalized[accession]
        output = pair_outputs[(gene, accession)]
        selected_hit = selection.get(gene, {}) if selection.get(gene, {}).get("accession") == accession else {}
        pair_source = selected_hit.get("source", "direct" if accession in direct_by_gene[gene] else "none")
        all_rows.extend({**_tsv_row(gene, accession, ligand),
                         "selection_source": pair_source,
                         "mapped_uniprot": ";".join(sorted(direct_by_gene[gene])),
                         "sequence_identity": selected_hit.get("identity", ""),
                         "query_coverage": selected_hit.get("query_coverage", ""),
                         "subject_coverage": selected_hit.get("subject_coverage", ""),
                         "cif_path": output["cif_path"], "pymol_script_path": output["pymol_script_path"]}
                        for ligand in ligands)
        accepted = [ligand for ligand in ligands if ligand["filter_status"] == "accepted"]
        status = statuses[accession]
        usable = status.get("metadata_usable", False)
        hit = selection.get(gene, {}) if selection.get(gene, {}).get("accession") == accession else {}
        source = hit.get("source", "direct" if accession in direct_by_gene[gene] else "none")
        summary_rows.append({
            "gene": gene, "uniprot": accession, "selection_source": source,
            "mapped_uniprot": ";".join(sorted(direct_by_gene[gene])),
            "sequence_identity": hit.get("identity", ""),
            "query_coverage": hit.get("query_coverage", ""),
            "subject_coverage": hit.get("subject_coverage", ""),
            "download_source": status.get("download_source", ""),
            "alphafill_status": status.get("status"),
            "metadata_usable": usable, "total_transfers": len(ligands),
            "accepted_transfers": len(accepted), "rejected_transfers": len(ligands) - len(accepted),
            **{f"{category}_confidence_transfers": sum(
                item["confidence_category"] == category for item in ligands
            ) for category in CONFIDENCE_CATEGORIES},
            **summarize_compounds(ligands, usable=usable), **output,
            "metadata_path": status.get("metadata_path", ""),
            "error_message": status.get("error_message", ""),
        })
        genes_json.setdefault(gene, {})[accession] = {
            "download": status, "model_output": output, "thresholds": thresholds, "ligands": ligands,
        }

    all_df = pd.DataFrame(all_rows, columns=LIGAND_COLUMNS).sort_values(["gene", "uniprot", "ligand_instance_id"]) if all_rows else pd.DataFrame(columns=LIGAND_COLUMNS)
    filtered_df = all_df[all_df["filter_status"] == "accepted"].copy()
    category_columns = ["gene", "uniprot", "confidence_category", "identity_category", "transfer_count"]
    if all_df.empty:
        category_counts = pd.DataFrame(columns=category_columns)
    else:
        category_counts = (all_df.groupby(category_columns[:-1], sort=True)
                           .size().reset_index(name="transfer_count"))
    accession_summary = pd.DataFrame(summary_rows, columns=SUMMARY_COLUMNS).sort_values(["gene", "uniprot"]) if summary_rows else pd.DataFrame(columns=SUMMARY_COLUMNS)
    summary_df = compact_summary(accession_summary)
    failures = _execution_failures(statuses, normalized, pair_outputs)
    for filename, frame in (("alphafill_ligands_all.tsv", all_df),
                            ("alphafill_ligands_filtered.tsv", filtered_df),
                            ("alphafill_summary.tsv", summary_df),
                            ("alphafill_accession_summary.tsv", accession_summary),
                            ("alphafill_category_counts.tsv", category_counts)):
        _atomic_bytes(root / filename, _portable_frame(frame, root).to_csv(sep="\t", index=False).encode())
    _atomic_json(root / "alphafill_results.json", _portable_paths({
        "path_base": ".",
        "path_description": "All artifact *_path fields are relative to the ligands/alphafill directory containing this report",
        "execution": {"status": "completed_with_errors" if failures else "completed",
                      "strict_errors": options["strict_errors"], "failure_count": len(failures),
                      "failures": failures},
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "thresholds": thresholds,
        "category_definitions": {
            "confidence": {"high": "local_rmsd < high_local_rmsd AND tcs < high_tcs",
                           "intermediate": "not high; local_rmsd < intermediate_local_rmsd AND tcs < intermediate_tcs",
                           "low": "valid metrics outside high/intermediate limits",
                           "unknown": "missing or invalid local_rmsd or tcs"},
            "identity": {"below_25": "0 <= identity < 25",
                         **{f"{lower}-{lower + 10}": f"{lower} <= identity < {lower + 10}"
                            for lower in range(25, 75, 10)},
                         ">=75": "75 <= identity <= 100", "unknown": "missing or invalid identity"},
            "global_rmsd": "reported only; never used for classification",
            "accepted": "high or intermediate confidence, independently of identity",
        },
        "compound_summary_rule": "Unique placed compound codes per locus across all UniProt accessions; best known confidence wins; high + intermediate = accepted unique count",
        "selection_rule": "Direct mapped accessions first; if none has transplants, use only the highest-ranked BLAST candidate meeting both sequence coverage limits, whether or not it has transplants",
        "blast_file": str(blast_path) if blast_path else "",
        "gene_summaries": json.loads(summary_df.to_json(orient="records")),
        "genes": genes_json,
    }, root))
    _atomic_bytes(root / "README.txt", (
        "AlphaFill results from FastTarget\n\n"
        "Start with alphafill_selection.tsv: one row per locus. mapped_uniprot is the "
        "original mapping; selected_uniprot is the chosen evidence source. "
        "selection_source is direct, exact (100% sequence match), similar, or none. "
        "A fallback uses only the top BLAST candidate with identity and coverage of "
        "both proteins at or above the configured thresholds; it stops even if that "
        "candidate has no transplants.\n\n"
        "result_status: transplants_found = AlphaFill entry with placements; "
        "entry_found_without_transplants = valid entry with zero placements; "
        "alphafill_entry_not_found = HTTP 404; incomplete_download = network, "
        "schema or file error; unmapped_or_no_candidates = no eligible UniProt ID.\n\n"
        "alphafill_summary.tsv has one row per locus for the FastTarget merge. "
        "alphafill_accession_summary.tsv has one row per locus and evaluated accession, "
        "with transfer counts, selection provenance, download source and file paths. "
        "alphafill_ligands_all.tsv lists every placed compound; "
        "alphafill_ligands_filtered.tsv lists accepted high/intermediate confidence "
        "placements. alphafill_category_counts.tsv counts placements by category. "
        "alphafill_results.json contains the same data and execution failures.\n\n"
        "Path columns are relative to this directory. Raw rsync files live in "
        "rsync_raw/. Valid downloads and analysis checkpoints live in cache/. "
        "Reruns reuse them; refresh_cache: True requests fresh data.\n"
    ).encode())
    if options["strict_errors"] and failures:
        raise AlphaFillExecutionError(failures, root / "alphafill_results.json")
    return summary_df


def compact_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Combine accession summaries into one row per locus, without double counts.

    Merge compound lists, never sum their counts: the same compound may occur
    in several UniProt accessions or confidence categories. A partial result
    reports observed compounds and alphafill_complete=False; unavailable data
    has null counts rather than implying a successful zero-hit search.
    """
    if summary.empty:
        return pd.DataFrame(columns=GENE_SUMMARY_COLUMNS)
    required = {"gene", "alphafill_status", "metadata_usable"} | {
        f"alphafill_{category}_confidence_compounds" for category in CONFIDENCE_CATEGORIES
    }
    missing = required - set(summary.columns)
    if missing:
        raise ValueError("Missing AlphaFill compound summary columns: " + ", ".join(sorted(missing)))
    rows = []
    for gene, group in summary.groupby("gene", sort=True):
        ligands = []
        for category in CONFIDENCE_CATEGORIES:
            for value in group[f"alphafill_{category}_confidence_compounds"]:
                if isinstance(value, str):
                    ligands.extend({"compound_code": code, "confidence_category": category}
                                   for code in value.split(";") if code.strip())
        rows.append({
            "gene": gene,
            "alphafill_status": ";".join(sorted(set(group["alphafill_status"].dropna().astype(str)))),
            "alphafill_complete": bool(group["metadata_usable"].fillna(False).all()),
            **summarize_compounds(ligands, usable=bool(group["metadata_usable"].fillna(False).any())),
        })
    result = pd.DataFrame(rows, columns=GENE_SUMMARY_COLUMNS)
    for column in COMPOUND_SUMMARY_COLUMNS:
        if column.endswith("_count") or column in {"alphafill_unique_compounds", "alphafill_unique_accepted_compounds"}:
            result[column] = result[column].astype("Int64")
    return result
