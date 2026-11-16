"""eggNOG-mapper annotations (COG, KEGG, Pfam, CAZy) keyed by genome locus tag.

The user runs eggNOG-mapper (e.g. the eggNOG 7 web server) on the organism .faa
and gives the ``<prefix>.emapper.annotations`` file; nothing is computed here.

Confidence (emapper >= 3): ``annotation_confidence`` holds one character per
functional column, h/m/l, from the orthology type of the donor the value was
copied from (h: seed or one-to-one ortholog, m: one-to-many or many-to-one,
l: many-to-many). Values below ``min_confidence`` are dropped and their field is
listed in ``eggnog_low_confidence_dropped``, so the absence is explained.

COG: eggNOG 7 writes the COG id (COG0347) in ``COG_category``; older versions
write the functional letters. Both are resolved with the NCBI COG 2024 tables.
"""

import io
import logging
from pathlib import Path

import pandas as pd

from ftscripts import metadata
from ftscripts.functional_annotations import GENERIC_MAPS

DATA_DIR = Path(__file__).resolve().parent / "data"
COG_DEF_FILE = DATA_DIR / "cog-24.def.tab"
COG_FUN_FILE = DATA_DIR / "cog-24.fun.tab"

# Documented default of emapper 3; the file header overrides it when present.
DEFAULT_CONFIDENCE_ORDER = ("Preferred_name GOs EC KEGG_ko KEGG_Pathway KEGG_Module KEGG_Reaction "
                            "KEGG_rclass BRITE KEGG_TC CAZy BiGG_Reaction PFAMs").split()
CONFIDENCE_CODES = {"h": "high", "m": "medium", "l": "low"}
CONFIDENCE_LEVELS = ("high", "medium", "low")
DEFAULT_MIN_CONFIDENCE = "medium"

# eggNOG column -> output column. GO is left out: the GO columns come from UniProt.
# EC is left out: it can differ from the EC of the metabolic model (e.g. EC numbers
# moved to class 7 since 2018) and be read as a model annotation. The model EC per
# gene and reaction is in the MetaGraphTools outputs (EC_Number column).
FIELDS = {
    "Preferred_name": "eggnog_preferred_name",
    "KEGG_ko": "eggnog_kegg_ko",
    "KEGG_Pathway": "eggnog_kegg_pathway",
    "KEGG_Module": "eggnog_kegg_module",
    "CAZy": "eggnog_cazy",
    "PFAMs": "eggnog_pfams",
}
COG_COLUMNS = ("eggnog_cog", "eggnog_cog_name", "eggnog_cog_category",
               "eggnog_cog_category_desc", "eggnog_cog_pathway")
OUTPUT_COLUMNS = (("gene", "eggnog_seed_ortholog", "eggnog_evalue", "eggnog_description")
                  + COG_COLUMNS + tuple(FIELDS.values())
                  + ("eggnog_medium_confidence_fields", "eggnog_low_confidence_dropped", "eggnog_og"))
# Fewer than this fraction of the file's queries among the genome loci means
# the file belongs to another organism or another annotation.
MIN_LOCUS_MATCH = 0.5


def read_annotations(path):
    """Return (table, confidence field order) from an .emapper.annotations file."""
    order = None
    body = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("## confidence field order:"):
                order = line.split(":", 1)[1].split()
            elif not line.startswith("##"):
                body.append(line)
    if not body:
        raise ValueError(f"{path} has no annotation table")
    table = pd.read_csv(io.StringIO("".join(body)), sep="\t", dtype=str, keep_default_na=False)
    table.columns = [column.lstrip("#") for column in table.columns]
    query = next((c for c in ("query", "query_name") if c in table.columns), None)
    if query is None:
        raise ValueError(f"{path} is not an eggNOG-mapper annotations file (no query column)")
    table = table.rename(columns={query: "gene"})
    return table, order or list(DEFAULT_CONFIDENCE_ORDER)


def load_cog_tables(def_file=COG_DEF_FILE, fun_file=COG_FUN_FILE):
    """COG id -> (letters, name, pathway) and letter -> category description."""
    cogs = {}
    with open(def_file, encoding="latin-1") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5 and parts[0].startswith("COG"):
                cogs[parts[0]] = (parts[1], parts[2], parts[4])
    letters = {}
    with open(fun_file, encoding="latin-1") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and len(parts[0]) == 1 and parts[0].isalpha():
                letters[parts[0]] = parts[-1]
    return cogs, letters


def _values(text):
    """Split an eggNOG multi-value cell; '-' and empty mean no value."""
    if not text or text == "-":
        return []
    return [value.strip() for value in text.split(",") if value.strip() and value.strip() != "-"]


def _joined(values):
    return "; ".join(dict.fromkeys(values))


def resolve_cog(text, cogs, letters):
    """COG columns for one gene from COG ids (eggNOG 7) or category letters (older)."""
    ids, category_letters, names, pathways = [], [], [], []
    for value in _values(text):
        if value.startswith("COG"):
            ids.append(value)
            category, name, pathway = cogs.get(value, ("", "", ""))
            category_letters.extend(category)
            names.append(name)
            pathways.append(pathway)
        elif value.isalpha():
            category_letters.extend(value)
    category_letters = list(dict.fromkeys(category_letters))
    return {
        "eggnog_cog": _joined(ids),
        "eggnog_cog_name": _joined(n for n in names if n),
        "eggnog_cog_category": "".join(category_letters),
        "eggnog_cog_category_desc": _joined(letters.get(c, c) for c in category_letters),
        "eggnog_cog_pathway": _joined(p for p in pathways if p),
    }


def _field_value(field, text):
    values = _values(text)
    if field == "KEGG_ko":
        values = [v.split(":", 1)[-1] for v in values]
    if field == "KEGG_Pathway":
        # Same convention as functional_annotations: drop the global overview maps,
        # and accept both bare numbers (eggNOG 7) and ko/map prefixes (older).
        numbers = [v[-5:] for v in values if v[-5:].isdigit()]
        values = [n for n in numbers if n not in GENERIC_MAPS]
    return _joined(values)


def annotate_row(row, order, min_confidence, cogs, letters):
    allowed = CONFIDENCE_LEVELS[:CONFIDENCE_LEVELS.index(min_confidence) + 1]
    codes = row.get("annotation_confidence", "")
    confidence = {field: CONFIDENCE_CODES.get(code) for field, code in zip(order, codes)} if codes else {}
    record = {
        "gene": row["gene"],
        "eggnog_seed_ortholog": row.get("seed_ortholog", ""),
        "eggnog_evalue": row.get("evalue", ""),
        "eggnog_description": "" if row.get("Description", "-") == "-" else row.get("Description", ""),
        **resolve_cog(row.get("COG_category", ""), cogs, letters),
        "eggnog_og": _joined(_values(row.get("eggNOG_OGs", ""))),
    }
    medium, dropped = [], []
    for field, column in FIELDS.items():
        value = _field_value(field, row.get(field, ""))
        level = confidence.get(field)
        if value and level is not None and level not in allowed:
            dropped.append(column.removeprefix("eggnog_"))
            value = ""
        elif value and level == "medium":
            medium.append(column.removeprefix("eggnog_"))
        record[column] = value
    record["eggnog_medium_confidence_fields"] = "; ".join(medium)
    record["eggnog_low_confidence_dropped"] = "; ".join(dropped)
    return record


def build_table(annotations_file, loci, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """One row per genome locus; loci without an eggNOG hit have empty columns."""
    if min_confidence not in CONFIDENCE_LEVELS:
        raise ValueError(f"min_confidence must be one of {CONFIDENCE_LEVELS}, got {min_confidence!r}")
    table, order = read_annotations(annotations_file)
    if "annotation_confidence" not in table.columns:
        logging.warning(f"{annotations_file} has no annotation_confidence column (eggNOG-mapper < 3): "
                        f"all values are kept without a confidence filter")
    queries = set(table["gene"])
    known = queries & set(loci)
    if queries and len(known) / len(queries) < MIN_LOCUS_MATCH:
        raise ValueError(f"Only {len(known)} of {len(queries)} queries in {annotations_file} are genome "
                         f"locus tags: was eggNOG-mapper run on this organism's .faa?")
    if len(known) < len(queries):
        extra = sorted(queries - known)
        logging.warning(f"{len(extra)} eggNOG queries are not genome locus tags and are ignored, "
                        f"e.g. {', '.join(extra[:5])}")
    cogs, letters = load_cog_tables()
    records = {row["gene"]: annotate_row(row, order, min_confidence, cogs, letters)
               for row in table.to_dict("records") if row["gene"] in known}
    empty = {column: "" for column in OUTPUT_COLUMNS}
    rows = [records.get(locus, {**empty, "gene": locus}) for locus in loci]
    table = pd.DataFrame(rows, columns=list(OUTPUT_COLUMNS))
    if not table["eggnog_description"].astype(bool).any():
        # Only eggNOG-mapper < 3 writes a Description column.
        table = table.drop(columns="eggnog_description")
    # Missing values as NA, so the run summary counts genes with a value correctly.
    return table.mask(table.eq(""))


def annotate(output_path, organism_name, annotations_file, min_confidence=DEFAULT_MIN_CONFIDENCE):
    """Write metabolism/eggnog_annotations.tsv and return it for FastTarget's merge."""
    loci = metadata.ref_gbk_locus(output_path, organism_name)
    table = build_table(annotations_file, loci, min_confidence)
    destination = Path(output_path) / organism_name / "metabolism" / "eggnog_annotations.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, sep="\t", index=False)
    annotated = table["eggnog_seed_ortholog"].notna().sum()
    dropped = table["eggnog_low_confidence_dropped"].notna().sum()
    logging.info(f"eggNOG annotations: {annotated} of {len(table)} genes with a hit; "
                 f"{dropped} genes had values below {min_confidence} confidence removed -> {destination}")
    return table
