"""Standalone local ChEMBL extraction, independent of pipeline integration.

Use the fasttarget conda environment. No REST API, DuckDB or extra packages.

    python -m ftscripts.chembl --database /path/to/chembl_37.db --taxid 287 \
        --taxonomy-dir /path/to/ncbi_taxonomy --proteins /path/to/proteins.faa \
        --output-dir /path/to/organism/ligands/chembl --include-human

Taxonomy expansion uses NCBI nodes.dmp (merged.dmp and names.dmp recommended).
An explicit --taxids-file is supported but never claimed to be complete.
Bacterial targets are SINGLE PROTEIN. Original and parental compound IDs are
preserved; summary counts deduplicate parental ChEMBL IDs. Activity categories
are not inferred. Missing structures and missing pChEMBL do not drop activities.

ChEMBL assignment (D/H/other) and our genome mapping are separate evidence axes.
The original ChEMBL homology identity is NOT supplied and remains unknown.
Mapping defaults: identity >=90%, coverage >=70% of BOTH donor and genome
protein, best passing bit score then E-value per donor, all exact ties retained.
Known variant sequences are used; missing variant sequences are explicitly
marked as reference-only, never presented as exact experimental sequence matches.
Human evidence includes proteins and complexes, without assigning complex
activity to subunits. Species context includes structured cellular/organism
assays and other molecular targets. Organism evidence is not attributed to a
protein; records without an explicit assay taxid cannot establish a cellular
experiment in the species. Outputs are TSV/CSV
and run_metadata.json; the locus summary has five columns, as for BindingDB.

Only chembl_ligands and the context detail tables contain the extra evidence.
Species assay counts deduplicate assay_id, not activity rows. Cellular and
whole-organism categories may overlap. species_molecular_target_count covers
all assigned molecular targets, including the original targets; the relative
other-target table excludes strain/variant donors mapping to the same locus.
Targets that cannot be mapped to the genome remain explicitly labeled candidates.
No activity/inactivity conclusion is derived from the presence of an assay.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import logging
from pathlib import Path
import re
import sqlite3
import subprocess
import sys

import pandas as pd
from Bio import SeqIO

# Pure shared utilities; importing these does not open or run BindingDB.
from ftscripts.bindingdb import expand_taxids, select_best_hits

LOG = logging.getLogger(__name__)
HUMAN_TYPES = ("SINGLE PROTEIN", "PROTEIN COMPLEX", "PROTEIN COMPLEX GROUP", "PROTEIN NUCLEIC-ACID COMPLEX")
SUMMARY_COLUMNS = ["gene", "chembl_status", "chembl_target_count", "chembl_compounds", "chembl_compound_count"]
DONOR_COLUMNS = ["donor_id", "tid", "target_chembl_id", "component_id", "accession", "tax_id", "variant_id",
                 "mutation", "sequence", "sequence_basis"]
MAPPING_COLUMNS = ["gene", "donor_id", "tid", "component_id", "variant_id", "chembl_taxid", "identity",
                   "genome_coverage", "chembl_coverage", "evalue", "bitscore", "exact_sequence_match",
                   "sequence_basis", "best_hit_tie_count"]
LEVEL_COLUMNS = ["gene", "evidence_level", "chembl_compounds", "chembl_compound_count", "measurement_count"]
LIGAND_COLUMNS = ["parent_molregno", "compound_id", "pref_name", "canonical_smiles", "standard_inchi_key",
                  "original_compound_ids", "bacterial_target_count", "bacterial_measurement_count",
                  "bacterial_direct_measurement_count", "bacterial_homology_measurement_count",
                  "human_target_count", "human_protein_target_count", "human_complex_target_count",
                  "human_measurement_count", "human_direct_measurement_count", "human_homology_measurement_count",
                  "species_cellular_assay_count", "species_cellular_no_molecular_target_assay_count",
                  "species_organism_assay_count", "species_no_molecular_target_assay_count",
                  "species_molecular_target_count", "species_molecular_target_ids"]
REQUIRED = {
    "target_dictionary": {"tid", "chembl_id", "tax_id", "target_type", "pref_name"},
    "target_components": {"tid", "component_id", "targcomp_id", "homologue"},
    "component_sequences": {"component_id", "accession", "sequence", "tax_id", "component_type"},
    "activities": {"activity_id", "assay_id", "molregno", "standard_type", "standard_value", "standard_relation",
                   "standard_units", "pchembl_value", "data_validity_comment", "potential_duplicate", "doc_id"},
    "assays": {"assay_id", "tid", "chembl_id", "description", "doc_id", "relationship_type", "confidence_score",
               "assay_tax_id", "assay_organism", "assay_strain", "assay_type", "variant_id", "bao_format", "cell_id", "assay_cell_type"},
    "target_type": {"target_type", "parent_type"}, "bioassay_ontology": {"bao_id", "label"},
    "molecule_dictionary": {"molregno", "chembl_id", "pref_name"},
    "molecule_hierarchy": {"molregno", "parent_molregno"},
    "compound_structures": {"molregno", "canonical_smiles", "standard_inchi_key"},
    "variant_sequences": {"variant_id", "mutation", "accession", "sequence"},
    "docs": {"doc_id", "chembl_id", "doi", "pubmed_id"}, "version": {"name"},
}


def _json(path, value):
    partial = path.with_name(path.name + ".part")
    partial.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    partial.replace(path)


def _export(root, name, frame):
    for extension, separator in (("tsv", "\t"), ("csv", ",")):
        path = root / f"{name}.{extension}"
        partial = path.with_name(path.name + ".part")
        frame.to_csv(partial, sep=separator, index=False)
        partial.replace(path)


def _frame(connection, sql, params=()):
    return pd.read_sql_query(sql, connection, params=params)


def _ids(connection, name, column, values):
    connection.execute(f"DROP TABLE IF EXISTS temp.{name}")
    connection.execute(f"CREATE TEMP TABLE {name} ({column} INTEGER PRIMARY KEY)")
    clean = sorted({int(v) for v in values if pd.notna(v)})
    connection.executemany(f"INSERT INTO {name} VALUES (?)", [(v,) for v in clean])


def inspect_schema(connection):
    schema = {}
    for table, required in REQUIRED.items():
        columns = [row[1] for row in connection.execute(f'PRAGMA table_info("{table}")')]
        if not required.issubset(columns):
            raise ValueError(f"Unsupported ChEMBL schema: {table} missing {sorted(required - set(columns))}")
        schema[table] = columns
    release = [row[0] for row in connection.execute("SELECT name FROM version WHERE name GLOB 'ChEMBL_[0-9]*'")]
    if not release:
        raise ValueError("No ChEMBL release identifier in the version table")
    return schema, release


def _targets(connection, taxids, target_types):
    _ids(connection, "target_taxids", "taxid", taxids)
    placeholders = ",".join("?" for _ in target_types)
    return _frame(connection, f"""SELECT t.* FROM target_dictionary t JOIN target_taxids s ON s.taxid=t.tax_id
                   WHERE t.target_type IN ({placeholders}) ORDER BY t.tid""", target_types)


def _components(connection, targets):
    _ids(connection, "component_targets", "tid", targets.tid.tolist())
    return _frame(connection, """SELECT tc.tid, tc.targcomp_id, tc.homologue, c.*
        FROM target_components tc JOIN component_targets q ON q.tid=tc.tid
        JOIN component_sequences c ON c.component_id=tc.component_id ORDER BY tc.tid,c.component_id""")


def _activities(connection, targets, parent_ids=None, *, species_taxids=None):
    """One row per activity; never join components onto measurements."""
    _ids(connection, "activity_targets", "tid", targets.tid.tolist())
    # Anchor the small selection with CROSS JOIN; otherwise SQLite can prefer
    # scanning the entire activities primary index to satisfy ORDER BY.
    source_join = "activity_targets q CROSS JOIN assays a ON a.tid=q.tid CROSS JOIN activities act ON act.assay_id=a.assay_id"
    if parent_ids is not None:
        _ids(connection, "selected_parents", "molregno", parent_ids)
        connection.execute("DROP TABLE IF EXISTS temp.selected_measured_molecules")
        connection.execute("CREATE TEMP TABLE selected_measured_molecules (molregno INTEGER PRIMARY KEY)")
        connection.execute("""INSERT OR IGNORE INTO selected_measured_molecules
            SELECT molregno FROM selected_parents UNION
            SELECT h.molregno FROM molecule_hierarchy h JOIN selected_parents p ON p.molregno=h.parent_molregno""")
        source_join = "selected_measured_molecules sm CROSS JOIN activities act ON act.molregno=sm.molregno CROSS JOIN assays a ON a.assay_id=act.assay_id JOIN activity_targets q ON q.tid=a.tid"
    where = ""
    if species_taxids is not None:
        if parent_ids is None:
            raise ValueError("Species context requires an explicit compound selection")
        _ids(connection, "context_taxids", "taxid", species_taxids)
        source_join = "selected_measured_molecules sm CROSS JOIN activities act ON act.molregno=sm.molregno CROSS JOIN assays a ON a.assay_id=act.assay_id"
        where = "WHERE a.assay_tax_id IN (SELECT taxid FROM context_taxids) OR t.tax_id IN (SELECT taxid FROM context_taxids)"
    result = _frame(connection, f"""SELECT act.*, a.tid, t.chembl_id AS target_chembl_id,
        t.target_type, t.pref_name AS target_name, t.tax_id AS target_taxid,
        a.chembl_id AS assay_chembl_id, a.description AS assay_description, a.doc_id AS assay_doc_id,
        a.assay_type, a.assay_organism, a.assay_tax_id, a.assay_strain, a.relationship_type,
        a.bao_format, b.label AS assay_format, a.cell_id, a.assay_cell_type,
        a.confidence_score, a.variant_id, v.mutation AS variant_mutation,
        v.accession AS variant_accession, v.sequence AS variant_sequence,
        md.chembl_id AS original_compound_id, h.parent_molregno AS reported_parent_molregno,
        coalesce(parent.molregno,md.molregno,act.molregno) AS parent_molregno,
        coalesce(parent.chembl_id,md.chembl_id) AS compound_id,
        coalesce(parent.pref_name,md.pref_name) AS compound_name,
        st.canonical_smiles, st.standard_inchi_key,
        CASE WHEN act.molregno IS NULL THEN 'missing_molecule_identifier'
             WHEN h.parent_molregno IS NOT NULL AND parent.molregno IS NULL THEN 'parent_dictionary_missing_used_original'
             WHEN md.molregno IS NULL THEN 'molecule_dictionary_missing'
             ELSE 'resolved' END AS compound_mapping_status
        FROM {source_join} LEFT JOIN target_dictionary t ON t.tid=a.tid
        LEFT JOIN bioassay_ontology b ON b.bao_id=a.bao_format
        LEFT JOIN molecule_dictionary md ON md.molregno=act.molregno
        LEFT JOIN molecule_hierarchy h ON h.molregno=act.molregno
        LEFT JOIN molecule_dictionary parent ON parent.molregno=h.parent_molregno
        LEFT JOIN compound_structures st ON st.molregno=coalesce(parent.molregno,md.molregno,act.molregno)
        LEFT JOIN variant_sequences v ON v.variant_id=a.variant_id {where} ORDER BY act.activity_id""")
    result["assignment_origin"] = result.relationship_type.map({"D": "direct", "H": "homology"}).fillna("other_or_unknown")
    result["chembl_homology_identity_pct"] = pd.Series(pd.NA, index=result.index, dtype="Float64")
    return result


def filter_activities(frame, *, direct_only=False, min_confidence_score=None,
                      require_pchembl=False, exclude_problematic=True):
    audit = frame.copy()
    reasons = []
    for row in audit.to_dict("records"):
        problems = []
        if direct_only and row["relationship_type"] != "D":
            problems.append("not_direct_assignment")
        score = row["confidence_score"]
        if min_confidence_score is not None and (pd.isna(score) or score < min_confidence_score):
            problems.append("confidence_below_minimum")
        if require_pchembl and pd.isna(row["pchembl_value"]):
            problems.append("missing_pchembl")
        validity = row["data_validity_comment"]
        if exclude_problematic and pd.notna(validity) and validity != "Manually validated":
            problems.append("problematic_data_validity")
        reasons.append(";".join(problems))
    audit["filter_reasons"] = reasons
    audit["filter_status"] = ["excluded" if value else "retained" for value in reasons]
    return audit, audit[audit.filter_status == "retained"].copy()


def species_context(connection, targets, parent_ids, taxids, *, direct_only=False,
                    min_confidence_score=None, require_pchembl=False, exclude_problematic=True):
    """Retrieve compound context using explicit assay/assigned-target taxonomy."""
    raw = _activities(connection, targets, parent_ids, species_taxids=taxids)
    hierarchy = dict(connection.execute("SELECT target_type,parent_type FROM target_type"))
    molecular = {"MOLECULAR"}
    while True:
        expanded = molecular | {kind for kind,parent in hierarchy.items() if parent in molecular}
        if expanded == molecular:
            break
        molecular = expanded
    raw["has_molecular_target"] = raw.target_type.isin(molecular)
    raw["assay_species_match"] = raw.assay_tax_id.isin(taxids)
    raw["assigned_target_species_match"] = raw.target_taxid.isin(taxids)
    raw["species_evidence_basis"] = [
        "assay_and_assigned_target_taxids" if a and t else "assay_taxid" if a else "assigned_target_taxid_only"
        for a,t in zip(raw.assay_species_match,raw.assigned_target_species_match)]
    raw["is_cellular_assay"] = (raw.assay_format.eq("cell-based format") | raw.cell_id.notna()
                                  | raw.assay_cell_type.fillna("").str.strip().ne(""))
    raw["is_whole_organism_assay"] = raw.assay_format.eq("organism-based format") | raw.target_type.eq("ORGANISM")
    raw["is_phenotypic_species_evidence"] = raw.assay_species_match & (raw.is_cellular_assay | raw.is_whole_organism_assay)
    raw["no_identified_molecular_target"] = ~raw.has_molecular_target
    # Assignment confidence and pChEMBL are protein-target filters, not valid
    # reasons to discard organism MIC or other phenotypic endpoints.
    molecular_audit,_ = filter_activities(raw[raw.has_molecular_target],direct_only=direct_only,
        min_confidence_score=min_confidence_score,require_pchembl=require_pchembl,exclude_problematic=exclude_problematic)
    nonmolecular_audit,_ = filter_activities(raw[~raw.has_molecular_target],exclude_problematic=exclude_problematic)
    audit = pd.concat([molecular_audit,nonmolecular_audit],ignore_index=True).sort_values("activity_id").reset_index(drop=True)
    retained = audit[audit.filter_status.eq("retained")].copy()
    phenotype = retained[retained.is_phenotypic_species_evidence].copy()
    phenotype["evidence_scope"] = "organism_or_cellular_not_attributed_to_protein"
    molecular_rows = retained[retained.has_molecular_target & retained.assigned_target_species_match].copy()
    return audit,retained,phenotype,molecular_rows


def other_target_context(pairs, mapping, molecular_rows, donors):
    """Other targets relative to each locus, not different strain entries of itself."""
    columns = ["gene","compound_id","other_tid","other_target_chembl_id","other_target_name",
               "other_target_type","other_mapped_genes","relationship_to_reference","assignment_origin",
               "species_evidence_basis","assay_species_match","assigned_target_species_match","measurement_count"]
    references = pairs[["gene","compound_id"]].dropna().drop_duplicates()
    related = _attach_donors(molecular_rows[molecular_rows.target_type.eq("SINGLE PROTEIN")],donors)
    remaining = molecular_rows[~molecular_rows.target_type.eq("SINGLE PROTEIN")].copy()
    remaining["donor_id"] = pd.array([pd.NA]*len(remaining),dtype="Int64")
    related = pd.concat([related,remaining],ignore_index=True)
    mapped = mapping.groupby("donor_id").gene.agg(set).to_dict()
    rows = []
    for reference in references.itertuples(index=False):
        for keys,group in related[related.compound_id.eq(reference.compound_id)].groupby(
                ["tid","assignment_origin","species_evidence_basis"],dropna=False,sort=True):
            tid,origin,basis = keys
            # A target may contain multiple assayed variants. Exclude only
            # the measurements whose donor maps back to the reference locus.
            keep = group.donor_id.map(lambda donor: reference.gene not in mapped.get(donor,set()) if pd.notna(donor) else True)
            group = group[keep]
            if group.empty:
                continue
            genes = set().union(*(mapped.get(d,set()) for d in group.donor_id.dropna()))
            first = group.iloc[0]
            unresolved = any(pd.isna(d) or not mapped.get(d,set()) for d in group.donor_id)
            relation = ("partially_mapped_target_candidate" if unresolved else "different_mapped_locus") if genes else (
                "unmapped_target_candidate" if first.target_type == "SINGLE PROTEIN" else "other_molecular_target_or_group")
            rows.append({"gene":reference.gene,"compound_id":reference.compound_id,"other_tid":int(tid),
                "other_target_chembl_id":first.target_chembl_id,"other_target_name":first.target_name,
                "other_target_type":first.target_type,"other_mapped_genes":";".join(sorted(genes)),
                "relationship_to_reference":relation,"assignment_origin":origin,"species_evidence_basis":basis,
                "assay_species_match":bool(group.assay_species_match.any()),"assigned_target_species_match":True,
                "measurement_count":group.activity_id.nunique()})
    return pd.DataFrame(rows,columns=columns)


def _sequence(value):
    return "" if pd.isna(value) else re.sub(r"\s+", "", str(value)).upper().rstrip("*")


def mapping_donors(targets, components, activities):
    rows = []
    for target in targets.itertuples():
        members = components[(components.tid == target.tid) & components.component_type.fillna("").str.upper().eq("PROTEIN")]
        if len(members) != 1:
            LOG.warning("Target %s has %s protein components; sequence assignment unresolved", target.chembl_id, len(members))
            continue
        component = members.iloc[0]
        combinations = [(None, None, None)]
        variants = activities[activities.tid == target.tid][["variant_id", "variant_mutation", "variant_sequence"]].drop_duplicates()
        for variant in variants.itertuples(index=False):
            if pd.notna(variant.variant_id):
                combinations.append((int(variant.variant_id), variant.variant_mutation, variant.variant_sequence))
        for variant_id, mutation, variant_sequence in combinations:
            variant_seq = _sequence(variant_sequence)
            basis = "variant_sequence" if variant_seq else ("canonical_sequence" if variant_id is None else "canonical_reference_variant_sequence_unavailable")
            rows.append({"donor_id": len(rows)+1, "tid": target.tid, "target_chembl_id": target.chembl_id,
                         "component_id": int(component.component_id), "accession": component.accession,
                         "tax_id": target.tax_id, "variant_id": variant_id, "mutation": mutation,
                         "sequence": variant_seq or _sequence(component.sequence), "sequence_basis": basis})
    return pd.DataFrame(rows, columns=DONOR_COLUMNS)


def map_proteins(proteins_fasta, donors, workdir, min_identity=90, min_coverage=70, threads=4):
    if not 0 < min_identity <= 100 or not 0 < min_coverage <= 100 or threads < 1:
        raise ValueError("Identity/coverage must be in (0,100]; threads must be positive")
    records = list(SeqIO.parse(str(proteins_fasta), "fasta"))
    genes = [record.id for record in records]
    if not records or len(genes) != len(set(genes)):
        raise ValueError("Protein FASTA must contain unique locus_tag identifiers")
    genome = {record.id: _sequence(str(record.seq)) for record in records}
    workdir.mkdir(parents=True, exist_ok=True)
    usable = donors[donors.sequence.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYBXZJUO]+")].copy()
    usable = usable[usable.sequence.map(lambda seq: bool(set(seq)-set("ACGTUN")))]
    query = workdir/"chembl_donors.faa"
    with query.open("w") as handle:
        for row in usable.itertuples():
            handle.write(f">CMBD{row.donor_id}\n{row.sequence}\n")
    if usable.empty:
        return pd.DataFrame(columns=MAPPING_COLUMNS), genes
    tools = {name: Path(sys.prefix)/"bin"/name for name in ("makeblastdb", "blastp")}
    for name, path in tools.items():
        if not path.is_file():
            raise FileNotFoundError(f"{name} missing from {sys.prefix}; add it to the fasttarget environment")
    prefix = workdir/"genome_proteins"
    subprocess.run([str(tools["makeblastdb"]), "-in", str(Path(proteins_fasta).resolve()), "-dbtype", "prot", "-out", str(prefix)], check=True, capture_output=True, text=True)
    raw = workdir/"blastp.tsv"
    subprocess.run([str(tools["blastp"]), "-query", str(query), "-db", str(prefix), "-out", str(raw),
                    "-outfmt", "6 qseqid sseqid pident qstart qend sstart send qlen slen evalue bitscore",
                    "-max_hsps", "1", "-max_target_seqs", str(max(5, len(genes))), "-evalue", "1e-5",
                    "-num_threads", str(threads)], check=True, capture_output=True, text=True)
    lookup = usable.set_index("donor_id")
    rows = []
    for line in raw.read_text().splitlines():
        q,s,identity,qs,qe,ss,se,qlen,slen,evalue,bitscore = line.split("\t")
        dcov = 100*(abs(int(qe)-int(qs))+1)/int(qlen)
        gcov = 100*(abs(int(se)-int(ss))+1)/int(slen)
        if float(identity) >= min_identity and min(dcov,gcov) >= min_coverage:
            donor_id = int(q.removeprefix("CMBD")); donor = lookup.loc[donor_id]
            rows.append({"gene":s, "donor_id":donor_id, "tid":int(donor.tid), "component_id":int(donor.component_id),
                         "variant_id":donor.variant_id, "chembl_taxid":int(donor.tax_id), "identity":float(identity),
                         "genome_coverage":gcov, "chembl_coverage":dcov, "evalue":float(evalue), "bitscore":float(bitscore),
                         "exact_sequence_match":donor.sequence == genome[s], "sequence_basis":donor.sequence_basis, "best_hit_tie_count":1})
    candidates = pd.DataFrame(rows, columns=MAPPING_COLUMNS)
    best,ties = select_best_hits(candidates.rename(columns={"donor_id":"polymerid"}), "all")
    best,ties = best.rename(columns={"polymerid":"donor_id"}), ties.rename(columns={"polymerid":"donor_id"})
    for frame in (best,ties):
        if not frame.empty:
            frame["best_hit_tie_count"] = frame.groupby("donor_id").gene.transform("size")
    candidates.to_csv(workdir/"chembl_passing_hits.tsv", sep="\t", index=False)
    ties.to_csv(workdir/"chembl_best_hit_ties.tsv", sep="\t", index=False)
    return best, genes


def _attach_donors(activities, donors):
    result = activities.copy()
    lookup = {(int(row.tid), int(row.variant_id) if pd.notna(row.variant_id) else None): row.donor_id for row in donors.itertuples()}
    result["donor_id"] = pd.array([lookup.get((int(row.tid), int(row.variant_id) if pd.notna(row.variant_id) else None)) for row in result.itertuples()], dtype="Int64")
    return result


def evidence_level(row):
    source = row["assignment_origin"]
    if row["sequence_basis"] == "canonical_reference_variant_sequence_unavailable":
        return source + "_variant_sequence_unresolved"
    return source + ("_identical_sequence" if row["exact_sequence_match"] else "_similar_sequence")


def _ligands(activities):
    rows = []
    for parent, group in activities[activities.compound_id.notna()].groupby("parent_molregno", sort=True):
        first = group.iloc[0]
        rows.append({"parent_molregno":int(parent), "compound_id":first.compound_id, "pref_name":first.compound_name,
                     "canonical_smiles":first.canonical_smiles, "standard_inchi_key":first.standard_inchi_key,
                     "original_compound_ids":";".join(sorted(set(group.original_compound_id.dropna()))),
                     "bacterial_target_count":group.tid.nunique(), "bacterial_measurement_count":group.activity_id.nunique(),
                     "bacterial_direct_measurement_count":group[group.assignment_origin.eq("direct")].activity_id.nunique(),
                     "bacterial_homology_measurement_count":group[group.assignment_origin.eq("homology")].activity_id.nunique()})
    return pd.DataFrame(rows, columns=LIGAND_COLUMNS)


def run_chembl(database_path, taxid, output_dir, *, proteins_fasta=None, taxonomy_dir=None, taxids=None,
               include_human=False, min_identity=90, min_coverage=70, threads=4, direct_only=False,
               min_confidence_score=None, require_pchembl=False, exclude_problematic=True, include_species_context=True):
    source = Path(database_path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if min_confidence_score is not None and not 0 <= min_confidence_score <= 9:
        raise ValueError("min_confidence_score must be between 0 and 9")
    expanded,taxonomy = expand_taxids(taxid, taxonomy_dir, taxids)
    root = Path(output_dir).expanduser().resolve(); root.mkdir(parents=True, exist_ok=True)
    metadata = {"status":"running", "taxid":int(taxid), "taxids":expanded, "taxonomy":taxonomy,
                "database":source.name, "started_at":datetime.now(timezone.utc).isoformat(),
                "filters":{"direct_only":direct_only, "min_confidence_score":min_confidence_score,
                           "require_pchembl":require_pchembl, "exclude_problematic":exclude_problematic},
                "mapping":{"min_identity":min_identity, "min_coverage_both_sequences":min_coverage,
                           "direction":"chembl_donor_to_genome", "selection":"highest_bitscore_then_lowest_evalue",
                           "ties":"assign_all_equal_best", "variants":"use_known_assay_sequence_else_explicit_reference_only"},
                "include_human":include_human, "human_target_types":list(HUMAN_TYPES), "organism_assays":"structured_species_context" if include_species_context else "not_searched",
                "species_context_policy":{"enabled":include_species_context,"phenotypic_taxonomy":"explicit_assay_taxid_required",
                    "phenotypic_filters":"data_validity_only; no protein confidence, direct-assignment or pChEMBL requirement",
                    "molecular_taxonomy":"assigned_target_taxid_preserved_separately_from_assay_taxid",
                    "other_target_identity":"exclude_donors_mapping_to_the_same_locus","no_positive_activity_inference":True}, "activity_categories":"not_assigned",
                "compound_identity":"parent_molregno with original fallback if parent dictionary missing",
                "chembl_homology_identity":"not_provided_by_source", "summary_scope":"union of retained evidence; levels in separate detail table"}
    metadata["database_file"] = {"bytes":source.stat().st_size, "mtime_ns":source.stat().st_mtime_ns}
    for directory in list(source.parents)[:4]:
        manifest = directory/"metadata.json"
        if manifest.is_file():
            upstream = json.loads(manifest.read_text())
            if upstream.get("source") == "ChEMBL" and (directory/upstream.get("database", "")).resolve() == source:
                metadata["download_provenance"] = {key:upstream.get(key) for key in
                    ("version", "source_url", "sha256", "checksum_verified")}
                break
    _json(root/"run_metadata.json", metadata)
    c = sqlite3.connect(source.as_uri()+"?mode=ro", uri=True)
    try:
        schema,release = inspect_schema(c); metadata["chembl_release"] = release
        _json(root/"source_schema.json", schema)
        targets = _targets(c, expanded, ("SINGLE PROTEIN",)); components = _components(c, targets)
        audit,activities = filter_activities(_activities(c, targets), direct_only=direct_only,
            min_confidence_score=min_confidence_score, require_pchembl=require_pchembl, exclude_problematic=exclude_problematic)
        donors = mapping_donors(targets, components, activities)
        activities = _attach_donors(activities, donors); ligands = _ligands(activities)
        LOG.info("Selected taxids=%s targets=%s activities=%s/%s ligands=%s", len(expanded),len(targets),len(activities),len(audit),len(ligands))
        if proteins_fasta:
            mapping,genes = map_proteins(proteins_fasta, donors, root/"mapping", min_identity, min_coverage, threads)
            metadata["mapping"]["proteins_sha256"] = hashlib.sha256(Path(proteins_fasta).read_bytes()).hexdigest()
        else:
            mapping,genes = pd.DataFrame(columns=MAPPING_COLUMNS), []
        pairs = mapping.merge(activities.drop(columns=["tid","variant_id"]), on="donor_id", how="inner")
        pairs["evidence_level"] = [evidence_level(row) for row in pairs.to_dict("records")]
        human_targets = targets.iloc[:0].copy(); human_components = components.iloc[:0].copy()
        human = audit.iloc[:0].copy(); human_audit = audit.iloc[:0].copy()
        cross = pd.DataFrame(columns=["compound_id","bacterial_tid","human_tid","human_target_chembl_id","human_target_type","human_assignment_origin"])
        if include_human:
            human_targets = _targets(c, [9606], HUMAN_TYPES)
            human_audit,human = filter_activities(_activities(c, human_targets, ligands.parent_molregno.tolist()),
                direct_only=direct_only, min_confidence_score=min_confidence_score, require_pchembl=require_pchembl,
                exclude_problematic=exclude_problematic)
            human_targets = human_targets[human_targets.tid.isin(human.tid)]; human_components = _components(c, human_targets)
            left = activities[["compound_id","tid"]].dropna(subset=["compound_id"]).drop_duplicates().rename(columns={"tid":"bacterial_tid"})
            right = human[["compound_id","tid","target_chembl_id","target_type","assignment_origin"]].drop_duplicates().rename(
                columns={"tid":"human_tid","target_chembl_id":"human_target_chembl_id","target_type":"human_target_type","assignment_origin":"human_assignment_origin"})
            cross = left.merge(right, on="compound_id", how="inner")
        for name, field, scope in (("human_target_count","tid",None), ("human_protein_target_count","tid","protein"),
                ("human_complex_target_count","tid","complex"), ("human_measurement_count","activity_id",None),
                ("human_direct_measurement_count","activity_id","direct"), ("human_homology_measurement_count","activity_id","homology")):
            subset = human
            if scope == "protein": subset = human[human.target_type.eq("SINGLE PROTEIN")]
            if scope == "complex": subset = human[human.target_type.ne("SINGLE PROTEIN")]
            if scope in {"direct","homology"}: subset = human[human.assignment_origin.eq(scope)]
            ligands[name] = ligands.parent_molregno.map(subset.groupby("parent_molregno")[field].nunique()).fillna(0).astype("Int64") if include_human else pd.Series(pd.NA,index=ligands.index,dtype="Int64")
        context_audit = audit.iloc[:0].copy()
        context = audit.iloc[:0].copy()
        phenotype = audit.iloc[:0].copy()
        molecular_context = audit.iloc[:0].copy()
        other_context = pd.DataFrame(columns=["gene","compound_id","other_tid","other_target_chembl_id",
            "other_target_name","other_target_type","other_mapped_genes","relationship_to_reference",
            "assignment_origin","species_evidence_basis","assay_species_match","assigned_target_species_match","measurement_count"])
        context_targets = targets.iloc[:0].copy(); context_components = components.iloc[:0].copy()
        if include_species_context:
            context_audit,context,phenotype,molecular_context = species_context(c,targets,ligands.parent_molregno.tolist(),expanded,
                direct_only=direct_only,min_confidence_score=min_confidence_score,require_pchembl=require_pchembl,exclude_problematic=exclude_problematic)
            other_context = other_target_context(pairs,mapping,molecular_context,donors)
            _ids(c,"context_target_ids","tid",context.tid.dropna())
            context_targets = _frame(c,"SELECT t.* FROM target_dictionary t JOIN context_target_ids s ON s.tid=t.tid")
            context_components = _components(c,context_targets)
        for name,frame in (("species_cellular_assay_count", phenotype[phenotype.is_cellular_assay] if include_species_context else phenotype),
                ("species_cellular_no_molecular_target_assay_count", phenotype[phenotype.is_cellular_assay & phenotype.no_identified_molecular_target] if include_species_context else phenotype),
                ("species_organism_assay_count", phenotype[phenotype.is_whole_organism_assay] if include_species_context else phenotype),
                ("species_no_molecular_target_assay_count", phenotype[phenotype.no_identified_molecular_target] if include_species_context else phenotype),
                ("species_molecular_target_count",molecular_context)):
            field = "tid" if name == "species_molecular_target_count" else "assay_id"
            ligands[name] = ligands.parent_molregno.map(frame.groupby("parent_molregno")[field].nunique()).fillna(0).astype("Int64") if include_species_context else pd.Series(pd.NA,index=ligands.index,dtype="Int64")
        target_lists = molecular_context.groupby("parent_molregno").target_chembl_id.agg(lambda values: ";".join(sorted(set(values.dropna()))))
        ligands["species_molecular_target_ids"] = ligands.parent_molregno.map(target_lists).fillna("") if include_species_context else pd.NA
        LOG.info("Species context: %s retained activities, %s cellular/organism rows, %s locus-compound-other-target rows",len(context),len(phenotype),len(other_context))
        summaries,levels = [],[]
        for gene in genes:
            mapped = mapping[mapping.gene.eq(gene)]; evidence = pairs[pairs.gene.eq(gene)]
            compounds = sorted(set(evidence.compound_id.dropna()))
            summaries.append({"gene":gene,"chembl_status":"no_target_mapping" if mapped.empty else ("no_activities" if evidence.empty else "found"),
                              "chembl_target_count":mapped.tid.nunique(),"chembl_compounds":";".join(compounds),"chembl_compound_count":len(compounds)})
            for level,group in evidence.groupby("evidence_level",sort=True):
                ids = sorted(set(group.compound_id.dropna()))
                levels.append({"gene":gene,"evidence_level":level,"chembl_compounds":";".join(ids),
                               "chembl_compound_count":len(ids),"measurement_count":group.activity_id.nunique()})
        summary = pd.DataFrame(summaries,columns=SUMMARY_COLUMNS)
        proteins = targets.rename(columns={"chembl_id":"target_chembl_id","tax_id":"target_taxid","organism":"target_organism"}).merge(components,on="tid",how="left")
        doc_ids,assay_ids = set(),set()
        for frame in (activities,human,context):
            doc_ids.update(frame.doc_id.dropna()); doc_ids.update(frame.assay_doc_id.dropna()); assay_ids.update(frame.assay_id)
        _ids(c,"reference_docs","doc_id",doc_ids); _ids(c,"selected_assays","assay_id",assay_ids)
        references = _frame(c,"SELECT d.* FROM docs d JOIN reference_docs r ON r.doc_id=d.doc_id")
        assays = _frame(c,"SELECT a.* FROM assays a JOIN selected_assays s ON s.assay_id=a.assay_id")
        variants = _frame(c,"SELECT DISTINCT v.* FROM variant_sequences v JOIN assays a ON a.variant_id=v.variant_id JOIN selected_assays s ON s.assay_id=a.assay_id")
        for name,frame in (("chembl_targets",targets),("chembl_proteins",proteins),("chembl_target_components",components),
                ("chembl_donors",donors),("chembl_activities_all",audit),("chembl_activities",activities),
                ("chembl_mapping",mapping),("chembl_protein_ligand_activities",pairs),("chembl_ligands",ligands),
                ("chembl_summary",summary),("chembl_evidence_levels",pd.DataFrame(levels,columns=LEVEL_COLUMNS)),
                ("chembl_assays",assays),("chembl_variants",variants),("chembl_references",references),
                ("chembl_human_targets",human_targets),("chembl_human_components",human_components),
                ("chembl_human_activities_all",human_audit),("chembl_human_activities",human),("chembl_human_crossactivity",cross),
                ("chembl_species_context_all",context_audit),("chembl_species_context",context),
                ("chembl_cellular_organism_activities",phenotype),("chembl_bacterial_molecular_activities",molecular_context),
                ("chembl_other_target_context",other_context),("chembl_context_targets",context_targets),("chembl_context_components",context_components)):
            _export(root,name,frame)
        metadata.update(status="completed",finished_at=datetime.now(timezone.utc).isoformat(),counts={
            "targets":len(targets),"donors":len(donors),"activities_before_filters":len(audit),"activities":len(activities),
            "direct_activities":int(activities.assignment_origin.eq("direct").sum()),"homology_activities":int(activities.assignment_origin.eq("homology").sum()),
            "ligands":len(ligands),"input_proteins":len(genes),"mapped_proteins":mapping.gene.nunique(),
            "human_activities":len(human),"human_complex_targets":human[human.target_type.ne("SINGLE PROTEIN")].tid.nunique(),
            "species_context_activities":len(context),"cellular_organism_activities":len(phenotype),"other_target_context_rows":len(other_context)})
        if targets.empty: LOG.warning("No SINGLE PROTEIN targets for the selected taxids")
        LOG.info("ChEMBL extraction complete: %s",metadata["counts"])
    except Exception as error:
        metadata.update(status="failed",error=str(error)); _json(root/"run_metadata.json",metadata)
        raise
    finally:
        c.close()
    _json(root/"run_metadata.json",metadata)
    return {"summary":summary,"metadata":metadata,"output_dir":root}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database",required=True); parser.add_argument("--taxid",type=int,required=True)
    parser.add_argument("--output-dir",required=True); parser.add_argument("--proteins")
    tax = parser.add_mutually_exclusive_group(); tax.add_argument("--taxonomy-dir"); tax.add_argument("--taxids-file")
    parser.add_argument("--include-human",action="store_true"); parser.add_argument("--direct-only",action="store_true")
    parser.add_argument("--min-confidence-score",type=int); parser.add_argument("--require-pchembl",action="store_true")
    parser.add_argument("--keep-problematic",action="store_true")
    parser.add_argument("--no-species-context",action="store_true",help="Skip compound cellular/organism/other-target context")
    parser.add_argument("--min-identity",type=float,default=90); parser.add_argument("--min-coverage",type=float,default=70)
    parser.add_argument("--threads",type=int,default=4)
    args = parser.parse_args(argv); logging.basicConfig(level=logging.INFO,format="%(asctime)s %(levelname)s: %(message)s")
    taxids = [int(v.strip()) for v in Path(args.taxids_file).read_text().splitlines() if v.strip()] if args.taxids_file else None
    run_chembl(args.database,args.taxid,args.output_dir,proteins_fasta=args.proteins,taxonomy_dir=args.taxonomy_dir,
        taxids=taxids,include_human=args.include_human,direct_only=args.direct_only,min_identity=args.min_identity,
        min_coverage=args.min_coverage,threads=args.threads,min_confidence_score=args.min_confidence_score,
        require_pchembl=args.require_pchembl,exclude_problematic=not args.keep_problematic,include_species_context=not args.no_species_context)


if __name__ == "__main__":
    main()
