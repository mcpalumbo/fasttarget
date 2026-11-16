"""Functional annotation joins must preserve locus tags and evidence."""

from ftscripts import functional_annotations as fa


SBML = '''<sbml xmlns:fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2"
 xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"><model>
 <fbc:listOfGeneProducts><fbc:geneProduct fbc:id="G_locusA" fbc:label="G_locusA">
  <rdf:li rdf:resource="https://identifiers.org/uniprot/P12345"/>
  <rdf:li rdf:resource="https://identifiers.org/kegg.genes/abc:otherA"/>
 </fbc:geneProduct></fbc:listOfGeneProducts>
 <listOfReactions><reaction id="R_one"><annotation>
  <rdf:li rdf:resource="https://identifiers.org/kegg.reaction/R00001"/>
 </annotation><fbc:geneProductAssociation><fbc:geneProductRef fbc:geneProduct="G_locusA"/>
 </fbc:geneProductAssociation></reaction></listOfReactions>
 </model></sbml>'''


def test_kegg_and_go_map_to_genome_locus_and_keep_evidence(tmp_path, monkeypatch):
    model = tmp_path / "model.xml"
    model.write_text(SBML)
    monkeypatch.setattr(fa.metadata, "ref_gbk_locus", lambda *_: ["locusA", "locusB"])
    responses = {
        "/get/br:br08901": "AMetabolism\nB  Carbohydrate metabolism\nC    00010  Glycolysis\nC    00020  TCA cycle\nC    01100  Metabolic pathways\n",
        "/list/pathway/abc": "path:abc00010\tGlycolysis\npath:abc00020\tTCA cycle\npath:abc01100\tMetabolic pathways\n",
        "/link/pathway/abc": "abc:otherA\tpath:abc00010\nabc:otherA\tpath:abc01100\n",
        "/link/pathway/rn": "rn:R00001\tpath:map00020\n",
        "organism_id%3A123": "Entry\tGene Names (ordered locus)\tGene Ontology (biological process)\tGene Ontology (molecular function)\tGene Ontology (cellular component)\nP12345\totherA\tmetabolic process [GO:0008152]\tcatalytic activity [GO:0003824]\t\n",
    }
    monkeypatch.setattr(fa, "fetch", lambda url: next(value for key, value in responses.items() if key in url))
    table = fa.annotate(tmp_path, "sample", 123, str(model), "abc", include_slim=False).set_index("gene")
    assert list(table.index) == ["locusA", "locusB"]
    assert "abc00010" in table.loc["locusA", "kegg_pathways"]
    assert "abc00020" in table.loc["locusA", "kegg_pathways"]
    assert "01100" not in table.loc["locusA", "kegg_pathways"]
    assert "reaction_inferred" in table.loc["locusA", "kegg_pathway_evidence"]
    assert "GO:0008152" in table.loc["locusA", "go_process"]
    assert table.loc["locusB", "kegg_pathways"] == ""


def test_go_slim_uses_ontology_ancestors():
    ontology = """[Term]
id: GO:0000001
name: broad
subset: goslim_generic

[Term]
id: GO:0000002
name: narrow
is_a: GO:0000001 ! broad
"""
    assert fa.go_slim_mapping(ontology)("GO:0000002") == {"broad [GO:0000001]"}


def test_kegg_requires_exact_taxon_match(monkeypatch):
    monkeypatch.setattr(fa, "fetch", lambda url: (
        "gn:T00001\ttaxid:10\n" if url.endswith("link/taxonomy/genome") else
        "T00001\tabc\tOrganism\n"))
    assert fa.kegg_organism(11) is None
    assert fa.kegg_organism(10) == "abc"
