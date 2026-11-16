"""Organism-independent KEGG and GO annotations keyed by genome locus tag."""

import csv
import io
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import pandas as pd

from ftscripts import metadata


KEGG = "https://rest.kegg.jp"
UNIPROT = "https://rest.uniprot.org/uniprotkb/stream"
GO_BASIC = "https://purl.obolibrary.org/obo/go/go-basic.obo"
RESOURCE = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"
GENERIC_MAPS = {"01100", "01110", "01120", "01200", "01210", "01212", "01220", "01230", "01232", "01240", "01250"}


def fetch(url):
    """Fetch a public annotation file; failures remain visible to the caller."""
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "FastTarget/functional-annotation"})
            with urllib.request.urlopen(request, timeout=90) as response:
                return response.read().decode("utf-8")
        except (OSError, TimeoutError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def pairs(text):
    for line in text.splitlines():
        if "\t" in line:
            yield line.split("\t", 1)


def joined(values):
    return "; ".join(sorted(set(values)))


def sbml_annotations(sbml_file):
    """Read explicit FBC gene associations and MIRIAM identifiers from any SBML."""
    genes = defaultdict(lambda: {"uniprot": set(), "kegg": set(), "reactions": set()})
    if not sbml_file:
        return genes
    root = ET.parse(sbml_file).getroot()
    product_ids = {}
    for element in root.iter():
        if element.tag.rsplit("}", 1)[-1] != "geneProduct":
            continue
        attributes = {key.rsplit("}", 1)[-1]: value for key, value in element.attrib.items()}
        product_id = attributes.get("id", "")
        locus = (attributes.get("label") or attributes.get("name") or product_id).removeprefix("G_")
        product_ids[product_id] = locus
        for node in element.iter():
            resource = node.get(RESOURCE, "")
            for source, field in (("uniprot/", "uniprot"), ("kegg.genes/", "kegg")):
                if source in resource:
                    genes[locus][field].add(urllib.parse.unquote(resource.rsplit(source, 1)[1]))
    for reaction in root.iter():
        if reaction.tag.rsplit("}", 1)[-1] != "reaction":
            continue
        reaction_id = reaction.get("id", "")
        for ref in reaction.iter():
            if ref.tag.rsplit("}", 1)[-1] == "geneProductRef":
                product = next((v for k, v in ref.attrib.items() if k.rsplit("}", 1)[-1] == "geneProduct"), "")
                locus = product_ids.get(product)
                if locus:
                    genes[locus]["reactions"].add(reaction_id)
    return genes


def kegg_organism(strain_taxid):
    """Resolve an exact NCBI strain taxid; never silently choose another strain."""
    genomes = {genome: taxon for genome, taxon in pairs(fetch(f"{KEGG}/link/taxonomy/genome"))}
    wanted = f"taxid:{strain_taxid}"
    matches = {genome.removeprefix("gn:") for genome, taxon in genomes.items() if taxon == wanted}
    codes = [line.split("\t")[1] for line in fetch(f"{KEGG}/list/organism").splitlines()
             if len(line.split("\t")) >= 2 and line.split("\t")[0] in matches]
    return codes[0] if len(codes) == 1 else None


def kegg_annotations(loci, sbml, strain_taxid, organism_code=None, sbml_file=None):
    """Assign organism-specific pathways and BRITE groups, retaining provenance."""
    code = organism_code or kegg_organism(strain_taxid)
    if not code:
        return {}, None
    hierarchy = {}
    major = minor = None
    for line in fetch(f"{KEGG}/get/br:br08901").splitlines():
        if line.startswith("A"):
            major = line[1:].strip()
        elif line.startswith("B"):
            minor = line[1:].strip()
        elif line.startswith("C"):
            match = re.search(r"\b(\d{5})\b", line)
            if match:
                hierarchy[match.group(1)] = (major, minor)
    names = {key[-5:]: name for key, name in pairs(fetch(f"{KEGG}/list/pathway/{code}"))}
    gene_paths = defaultdict(set)
    direct = defaultdict(set)
    for gene_id, pathway_id in pairs(fetch(f"{KEGG}/link/pathway/{code}")):
        locus = gene_id.split(":", 1)[-1]
        gene_paths[gene_id].add(pathway_id[-5:])
        if locus in loci:
            direct[locus].add(pathway_id[-5:])
    for locus, data in sbml.items():
        if locus not in loci:
            continue
        for gene_id in data["kegg"]:
            if gene_id.startswith(f"{code}:"):
                direct[locus].update(gene_paths[gene_id])
    # A reaction can occur in several pathways. Record these as inferred membership.
    reactions = defaultdict(set)
    if sbml_file and any(data["reactions"] for data in sbml.values()):
        root = ET.parse(sbml_file).getroot()
        if root is not None:
            for reaction in root.iter():
                if reaction.tag.rsplit("}", 1)[-1] != "reaction":
                    continue
                for node in reaction.iter():
                    resource = node.get(RESOURCE, "")
                    if "kegg.reaction/" in resource:
                        reactions[reaction.get("id", "")].add(resource.rsplit("kegg.reaction/", 1)[1])
                notes = " ".join(node.text or "" for node in reaction.iter()
                                 if node.tag.rsplit("}", 1)[-1] == "p")
                reactions[reaction.get("id", "")].update(re.findall(r"KEGG Reaction:\s*(R\d{5})", notes))
    reaction_paths = defaultdict(set)
    if reactions:
        for reaction_id, pathway_id in pairs(fetch(f"{KEGG}/link/pathway/rn")):
            reaction_paths[reaction_id.split(":", 1)[-1]].add(pathway_id[-5:])
    result = {}
    for locus in loci:
        inferred = {path for reaction in sbml.get(locus, {}).get("reactions", ())
                    for rxn in reactions[reaction] for path in reaction_paths[rxn] if path in names}
        specific = {path for path in direct[locus] | inferred
                    if path in hierarchy and hierarchy[path][0] == "Metabolism" and path not in GENERIC_MAPS}
        result[locus] = {
            "kegg_pathways": joined(f"{code}{path}: {names[path]}" for path in specific if path in names),
            "kegg_pathway_evidence": joined(f"{code}{path}: {'gene' if path in direct[locus] else 'reaction_inferred'}"
                                            for path in specific),
            "kegg_pathway_categories": joined(hierarchy[path][1] for path in specific if hierarchy[path][1]),
        }
    return result, code


def go_slim_mapping(obo_text):
    """Map GO terms to generic GO slim ancestors using is_a and part_of."""
    parents, slim, names = defaultdict(set), set(), {}
    term = None
    for line in obo_text.splitlines():
        if line == "[Term]":
            term = None
        elif line.startswith("id: GO:"):
            term = line[4:].strip()
        elif term and line.startswith("name: "):
            names[term] = line[6:].strip()
        elif term and line.startswith("is_a: GO:"):
            parents[term].add(line[6:16])
        elif term and line.startswith("relationship: part_of GO:"):
            parents[term].add(line[22:32])
        elif term and line == "subset: goslim_generic":
            slim.add(term)
    def ancestors(term):
        seen, pending = set(), [term]
        while pending:
            current = pending.pop()
            if current in seen:
                continue
            seen.add(current)
            pending.extend(parents[current] - seen)
        return {f"{names.get(value, value)} [{value}]" for value in seen & slim}
    return ancestors


def go_annotations(loci, sbml, strain_taxid, obo_text=None):
    """Use ordered locus tags and SBML UniProt accessions from the exact strain."""
    query = urllib.parse.urlencode({"query": f"organism_id:{strain_taxid}", "format": "tsv",
                                    "fields": "accession,gene_oln,go_p,go_f,go_c"})
    reader = csv.DictReader(io.StringIO(fetch(f"{UNIPROT}?{query}")), delimiter="\t")
    by_accession = {}
    by_locus = defaultdict(set)
    for row in reader:
        accession = row["Entry"]
        by_accession[accession] = row
        for locus in re.split(r"[\s;,]+", row.get("Gene Names (ordered locus)", "")):
            if locus in loci:
                by_locus[locus].add(accession)
    slim = go_slim_mapping(obo_text) if obo_text else None
    result = {}
    for locus in loci:
        accessions = by_locus[locus] | (sbml.get(locus, {}).get("uniprot", set()) & by_accession.keys())
        record = {"go_uniprot": joined(accessions)}
        for field, column, output in (("go_p", "Gene Ontology (biological process)", "go_process"),
                                      ("go_f", "Gene Ontology (molecular function)", "go_function"),
                                      ("go_c", "Gene Ontology (cellular component)", "go_component")):
            values = {term.strip() for accession in accessions for term in by_accession[accession][column].split(";")
                      if term.strip()}
            record[output] = joined(values)
            if slim:
                record[output + "_slim"] = joined(parent for term in values
                                                   for go_id in re.findall(r"GO:\d{7}", term)
                                                   for parent in slim(go_id))
        result[locus] = record
    return result


def annotate(output_path, organism_name, strain_taxid, sbml_file=None, organism_code=None, include_slim=True):
    """Write a one-row-per-genome-gene table and return it for FastTarget's merge."""
    loci = metadata.ref_gbk_locus(output_path, organism_name)
    sbml = sbml_annotations(sbml_file)
    kegg, code = kegg_annotations(set(loci), sbml, strain_taxid, organism_code, sbml_file)
    obo = fetch(GO_BASIC) if include_slim else None
    go = go_annotations(set(loci), sbml, strain_taxid, obo)
    rows = []
    for locus in loci:
        rows.append({"gene": locus, "kegg_organism": code or "",
                     "metabolic_model_gene": locus in sbml if sbml_file else None,
                     **kegg.get(locus, {"kegg_pathways": "", "kegg_pathway_evidence": "",
                                         "kegg_pathway_categories": ""}), **go[locus]})
    table = pd.DataFrame(rows)
    destination = Path(output_path) / organism_name / "metabolism" / "functional_annotations.tsv"
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, sep="\t", index=False)
    return table
