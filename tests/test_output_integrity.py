"""Outputs that exist but are not what the step wrote.

The parameters in a manifest say how an output was meant to be produced. They
cannot tell a complete file from one whose write was interrupted, so a table cut
short by a full disk or a killed process was reused as a valid cache entry.
"""

import io
import os
import contextlib

import pytest

from conftest import GENES, write_blast_output
from ftscripts import essentiality, files, genome, pathways


@pytest.fixture
def deg_table(organism):
    """A DEG table built from five hits, with its manifest."""
    output_path, organism_name = organism
    write_blast_output(
        os.path.join(output_path, organism_name, 'essentiality', 'deg_blast.tsv'),
        [(gene, 'DEG_1', 90.0, 95) for gene in GENES],
    )
    with contextlib.redirect_stdout(io.StringIO()):
        essentiality.deg_parse(output_path, organism_name, 50, 70)

    return output_path, organism_name, os.path.join(
        output_path, organism_name, 'essentiality', 'hit_in_deg.tsv'
    )


def test_a_complete_output_is_reused(deg_table):
    output_path, organism_name, path = deg_table

    with contextlib.redirect_stdout(io.StringIO()) as output:
        table = essentiality.deg_parse(output_path, organism_name, 50, 70)

    assert 'Rebuilding' not in output.getvalue()
    assert len(table) == len(GENES)


def test_a_truncated_output_is_rebuilt(deg_table):
    """Twenty rows cut to two used to be a cache hit."""
    output_path, organism_name, path = deg_table
    with open(path) as handle:
        lines = handle.readlines()
    with open(path, 'w') as handle:
        handle.writelines(lines[:3])

    with contextlib.redirect_stdout(io.StringIO()) as output:
        table = essentiality.deg_parse(output_path, organism_name, 50, 70)

    assert 'checksum' in output.getvalue()
    assert len(table) == len(GENES), 'the rebuilt table is complete again'


def test_an_edited_output_is_rebuilt(deg_table):
    """Same length, different content: the size check alone would miss this."""
    output_path, organism_name, path = deg_table
    with open(path) as handle:
        content = handle.read()
    with open(path, 'w') as handle:
        handle.write(content.replace('True', 'Fals'))

    assert files.manifest_output_changed(path) is True

    with contextlib.redirect_stdout(io.StringIO()) as output:
        essentiality.deg_parse(output_path, organism_name, 50, 70)

    assert 'checksum' in output.getvalue()


def test_an_output_with_no_recorded_checksum_is_unverifiable_not_broken(deg_table):
    """
    A manifest from an earlier version carries no checksum. Unknown is not the
    same as wrong, so the output is still reused.
    """
    output_path, organism_name, path = deg_table
    files.write_manifest(path, {'deg_identity_filter': 50, 'deg_coverage_filter': 70},
                         hash_output=False)

    assert files.manifest_output_changed(path) is False

    with contextlib.redirect_stdout(io.StringIO()) as output:
        essentiality.deg_parse(output_path, organism_name, 50, 70)

    assert 'Rebuilding' not in output.getvalue()


def test_the_checksum_and_size_are_recorded(deg_table):
    output_path, organism_name, path = deg_table

    manifest = files.read_manifest(path)

    assert manifest['output_sha256'] == files.sha256_file(path)
    assert manifest['output_bytes'] == os.path.getsize(path)


def test_a_missing_output_reports_as_changed(tmp_path):
    path = os.path.join(str(tmp_path), 'gone.tsv')
    with open(path, 'w') as handle:
        handle.write('content\n')
    files.write_manifest(path, {'a': 1})
    os.remove(path)

    assert files.manifest_output_changed(path) is True


def test_an_expensive_output_that_changed_refuses_rather_than_rebuilding(organism):
    """
    Roary and CoreCruncher are never re-run automatically, so a corrupt output
    stops the run and names the directory to remove.
    """
    from conftest import write_roary_inputs

    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)
    genes_csv = os.path.join(output_path, organism_name, 'conservation',
                             'roary_output', 'results', 'gene_presence_absence.csv')
    gff_dir = os.path.join(output_path, organism_name, 'conservation', 'roary_output', 'gff')
    files.write_manifest(genes_csv, {
        'min_core_freq': 99, 'min_identity': 95,
        'genome_count': len(os.listdir(gff_dir)),
        'genome_set': files.fileset_fingerprint(
            [os.path.join(gff_dir, name) for name in os.listdir(gff_dir)]
        ),
    })

    with open(genes_csv, 'a') as handle:
        handle.write('corrupted trailing line\n')

    with pytest.raises(RuntimeError, match='checksum'):
        genome.core_genome_programs(output_path, organism_name, core_threshold=99,
                                    identity=95, cpus=1, program_list=['roary'])
    assert os.path.exists(genes_csv), 'nothing is deleted by the refusal'


# --------------------------------------------------- degenerate metabolic networks

def _sif(output_path, organism_name, text):
    path = os.path.join(output_path, organism_name, 'metabolism', 'network.sif')
    with open(path, 'w') as handle:
        handle.write(text)
    return path


@pytest.mark.parametrize('content', ['', 'not a sif line\n', '   \n'])
def test_a_network_with_no_interactions_is_refused_by_name(organism, content):
    """
    max() over no components raised ValueError, which reached the caller as an
    unexplained failure of the whole metabolic module.
    """
    output_path, organism_name = organism
    path = _sif(output_path, organism_name, content)

    with pytest.raises(ValueError) as error:
        pathways.read_sif(output_path, organism_name)
    assert 'network.sif' in str(error.value)


@pytest.mark.parametrize('content,nodes', [
    ('A interacts B\n', 2),
    ('A interacts A\n', 1),
    ('A i B\nC i D\n', 2),
])
def test_a_component_too_small_to_have_a_centre_gives_zeros(organism, content, nodes):
    """
    Every centrality is zero when no node lies between two others, which is a
    legitimate result; dividing by that maximum raised ZeroDivisionError and lost
    the module.
    """
    output_path, organism_name = organism
    _sif(output_path, organism_name, content)

    with contextlib.redirect_stdout(io.StringIO()):
        graph, subgraph, centrality, degrees = pathways.read_sif(output_path, organism_name)

    assert subgraph.number_of_nodes() == nodes
    assert set(centrality.values()) == {0.0}
    assert degrees


def test_a_network_with_a_centre_is_normalised_to_one(organism):
    output_path, organism_name = organism
    _sif(output_path, organism_name, 'A i B\nB i C\nC i D\n')

    with contextlib.redirect_stdout(io.StringIO()):
        _, _, centrality, _ = pathways.read_sif(output_path, organism_name)

    assert max(centrality.values()) == 1.0
    assert min(centrality.values()) == 0.0


def test_only_the_largest_component_is_analysed(organism):
    output_path, organism_name = organism
    _sif(output_path, organism_name, 'A i B\nB i C\nC i D\nX i Y\n')

    with contextlib.redirect_stdout(io.StringIO()):
        graph, subgraph, centrality, _ = pathways.read_sif(output_path, organism_name)

    assert graph.number_of_nodes() == 6
    assert subgraph.number_of_nodes() == 4
    assert 'X' not in centrality
