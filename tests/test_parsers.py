"""Parsers on minimal inputs, including the boundaries.

The twenty-protein case exercises one point in each parameter space. These tests
pin the edges: a value exactly on the cutoff, no hits at all, a single genome,
and a threshold change over a cached table.
"""

import os

import pandas as pd
import pytest

from conftest import GENES, write_blast_output, write_corecruncher_inputs, write_roary_inputs
from ftscripts import essentiality, files, genome


def _deg_blast(output_path, organism_name, rows):
    write_blast_output(
        os.path.join(output_path, organism_name, 'essentiality', 'deg_blast.tsv'), rows
    )


def _true_genes(table, column):
    return sorted(table.loc[table[column].astype(str).str.lower() == 'true', 'gene'])


# --------------------------------------------------------------------------- DEG

def test_deg_keeps_only_hits_above_both_cutoffs(organism):
    output_path, organism_name = organism
    _deg_blast(output_path, organism_name, [
        (GENES[0], 'DEG_1', 80.0, 90),   # passes both
        (GENES[1], 'DEG_2', 80.0, 10),   # coverage too low
        (GENES[2], 'DEG_3', 20.0, 90),   # identity too low
    ])

    table = essentiality.deg_parse(output_path, organism_name, 50, 70)
    assert _true_genes(table, 'hit_in_deg') == [GENES[0]]


def test_deg_cutoffs_are_strict_so_a_value_on_the_boundary_is_excluded(organism):
    """The filter is '>' and not '>=', which only a boundary case reveals."""
    output_path, organism_name = organism
    _deg_blast(output_path, organism_name, [
        (GENES[0], 'DEG_1', 50.0, 71),   # identity exactly at the cutoff
        (GENES[1], 'DEG_2', 51.0, 70),   # coverage exactly at the cutoff
        (GENES[2], 'DEG_3', 51.0, 71),   # above both
    ])

    table = essentiality.deg_parse(output_path, organism_name, 50, 70)
    assert _true_genes(table, 'hit_in_deg') == [GENES[2]]


def test_deg_with_no_passing_hit_reports_every_gene_false(organism):
    output_path, organism_name = organism
    _deg_blast(output_path, organism_name, [(GENES[0], 'DEG_1', 10.0, 10)])

    table = essentiality.deg_parse(output_path, organism_name, 50, 70)
    assert _true_genes(table, 'hit_in_deg') == []
    assert len(table) == len(GENES), 'every gene must still have a row'


def test_deg_table_is_rebuilt_when_a_filter_changes(organism):
    output_path, organism_name = organism
    _deg_blast(output_path, organism_name, [
        (GENES[0], 'DEG_1', 95.0, 90),
        (GENES[1], 'DEG_2', 60.0, 90),
    ])

    strict = essentiality.deg_parse(output_path, organism_name, 90, 70)
    assert _true_genes(strict, 'hit_in_deg') == [GENES[0]]

    # Same output file, looser identity: the cached table would be wrong.
    relaxed = essentiality.deg_parse(output_path, organism_name, 50, 70)
    assert _true_genes(relaxed, 'hit_in_deg') == [GENES[0], GENES[1]]

    back = essentiality.deg_parse(output_path, organism_name, 90, 70)
    assert _true_genes(back, 'hit_in_deg') == [GENES[0]]


def test_deg_missing_blast_input_raises_naming_the_file(organism):
    output_path, organism_name = organism
    with pytest.raises(FileNotFoundError) as error:
        essentiality.deg_parse(output_path, organism_name, 50, 70)
    assert 'deg_blast.tsv' in str(error.value)


# ------------------------------------------------------------------------- Roary

def test_roary_core_is_the_genes_above_the_frequency_threshold(organism):
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {
        GENES[0]: 10, GENES[1]: 10, GENES[2]: 9, GENES[3]: 5, GENES[4]: 1,
    }, genomes=10)

    table = genome.roary_output(output_path, organism_name, core_threshold=0.99)
    assert _true_genes(table, 'core_roary') == [GENES[0], GENES[1]]


def test_roary_threshold_boundary_is_inclusive(organism):
    """The comparison is '>=', so a gene exactly at the threshold is core."""
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {
        GENES[0]: 10, GENES[1]: 9, GENES[2]: 8,
    }, genomes=10)

    table = genome.roary_output(output_path, organism_name, core_threshold=0.9)
    assert _true_genes(table, 'core_roary') == [GENES[0], GENES[1]]


def test_roary_table_is_rebuilt_when_the_threshold_changes(organism):
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {
        GENES[0]: 10, GENES[1]: 6, GENES[2]: 2,
    }, genomes=10)

    strict = genome.roary_output(output_path, organism_name, core_threshold=0.99)
    assert _true_genes(strict, 'core_roary') == [GENES[0]]

    relaxed = genome.roary_output(output_path, organism_name, core_threshold=0.5)
    assert _true_genes(relaxed, 'core_roary') == [GENES[0], GENES[1]]


def test_roary_single_genome_makes_every_present_gene_core(organism):
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 1, GENES[1]: 1}, genomes=1)

    table = genome.roary_output(output_path, organism_name, core_threshold=0.99)
    assert _true_genes(table, 'core_roary') == [GENES[0], GENES[1]]


def test_roary_missing_results_raises_instead_of_returning_an_unbound_name(organism):
    output_path, organism_name = organism
    with pytest.raises(FileNotFoundError):
        genome.roary_output(output_path, organism_name, core_threshold=0.99)


def test_a_cohort_that_does_not_match_the_prepared_genomes_is_refused(organism):
    """
    The denominator comes from the sample columns, so they have to be the genomes
    that were actually prepared: a GFF missing from the cohort would otherwise
    change the core genome without changing the reported threshold.
    """
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)
    gff_dir = os.path.join(output_path, organism_name, 'conservation', 'roary_output', 'gff')

    os.remove(os.path.join(gff_dir, 'genome_05.gff'))

    with pytest.raises(ValueError, match='prepared GFF cohort'):
        genome.roary_output(output_path, organism_name, core_threshold=0.99)


def test_a_cohort_without_the_reference_is_refused(organism):
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)
    gff_dir = os.path.join(output_path, organism_name, 'conservation', 'roary_output', 'gff')

    os.remove(os.path.join(gff_dir, f'{organism_name}.gff'))

    with pytest.raises(ValueError, match='prepared GFF cohort'):
        genome.roary_output(output_path, organism_name, core_threshold=0.99)


def test_an_unrecognised_presence_absence_schema_is_refused(organism):
    """A changed Roary output format must stop the parse, not be guessed at."""
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)
    csv_path = os.path.join(output_path, organism_name, 'conservation',
                            'roary_output', 'results', 'gene_presence_absence.csv')
    with open(csv_path, 'w') as handle:
        handle.write(f'"No. isolates","{organism_name}"\n10,"{GENES[0]}"\n')

    with pytest.raises(ValueError, match='schema'):
        genome.roary_output(output_path, organism_name, core_threshold=0.99)


def test_isolate_counts_must_agree_with_the_presence_cells(organism):
    """
    'No. isolates' is the number the threshold is applied to, so a value that
    disagrees with the row's own cells cannot be trusted.
    """
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 3}, genomes=10)
    csv_path = os.path.join(output_path, organism_name, 'conservation',
                            'roary_output', 'results', 'gene_presence_absence.csv')

    import csv as csv_module
    with open(csv_path, newline='') as handle:
        rows = list(csv_module.reader(handle))
    rows[1][3] = '10'          # claim ten isolates while only three cells are filled
    with open(csv_path, 'w', newline='') as handle:
        csv_module.writer(handle).writerows(rows)

    with pytest.raises(ValueError, match='disagree'):
        genome.roary_output(output_path, organism_name, core_threshold=0.99)


# ------------------------------------------------------------------ CoreCruncher

def test_corecruncher_core_comes_from_the_reference_entries(organism):
    output_path, organism_name = organism
    write_corecruncher_inputs(output_path, organism_name, [GENES[0], GENES[2]])

    table = genome.corecruncher_output(output_path, organism_name, core_threshold=99)
    assert _true_genes(table, 'core_corecruncher') == [GENES[0], GENES[2]]


def test_corecruncher_missing_input_raises_instead_of_declaring_no_core(organism):
    """
    A missing families_core.txt used to produce a table saying nothing is core,
    which is wrong data rather than a failure.
    """
    output_path, organism_name = organism
    os.makedirs(os.path.join(
        output_path, organism_name, 'conservation', 'corecruncher_output'
    ), exist_ok=True)

    with pytest.raises(FileNotFoundError) as error:
        genome.corecruncher_output(output_path, organism_name, core_threshold=99)
    assert 'families_core.txt' in str(error.value)


def test_corecruncher_empty_core_is_a_valid_result(organism):
    output_path, organism_name = organism
    write_corecruncher_inputs(output_path, organism_name, [])
    families = os.path.join(
        output_path, organism_name, 'conservation', 'corecruncher_output', 'families_core.txt'
    )
    with open(families, 'w') as handle:
        handle.write('other.faa&X\tother2.faa&Y\n')

    table = genome.corecruncher_output(output_path, organism_name, core_threshold=99)
    assert _true_genes(table, 'core_corecruncher') == []
    assert len(table) == len(GENES)


# ------------------------------------------------------- expensive-run refusals

def test_changed_core_parameters_refuse_to_reuse_an_expensive_run(organism):
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)
    write_corecruncher_inputs(output_path, organism_name, [GENES[0]])

    genes_csv = os.path.join(
        output_path, organism_name, 'conservation', 'roary_output', 'results',
        'gene_presence_absence.csv'
    )
    gff_dir = os.path.join(output_path, organism_name, 'conservation', 'roary_output', 'gff')
    recorded = {
        'min_core_freq': 99, 'min_identity': 95,
        'genome_count': len(os.listdir(gff_dir)),
        'genome_set': files.fileset_fingerprint(
            [os.path.join(gff_dir, name) for name in os.listdir(gff_dir)]
        ),
    }
    files.write_manifest(genes_csv, recorded)

    with pytest.raises(RuntimeError) as error:
        genome.core_genome_programs(output_path, organism_name, core_threshold=99, identity=90,
                                    cpus=1, program_list=['roary'])
    assert 'min_identity 95 -> 90' in str(error.value)
    assert os.path.exists(genes_csv), 'an expensive output must never be deleted'


def test_a_changed_genome_set_refuses_and_names_the_accession_file(organism):
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)

    genes_csv = os.path.join(
        output_path, organism_name, 'conservation', 'roary_output', 'results',
        'gene_presence_absence.csv'
    )
    files.write_manifest(genes_csv, {
        'min_core_freq': 99, 'min_identity': 95,
        'genome_count': 7, 'genome_set': 'deadbeefcafe',
    })

    with pytest.raises(RuntimeError) as error:
        genome.core_genome_programs(output_path, organism_name, core_threshold=99, identity=95,
                                    cpus=1, program_list=['roary'])
    message = str(error.value)
    assert 'genome_set' in message
    assert 'accession_file' in message


def test_an_unrecorded_parameter_warns_and_reuses(organism, capsys):
    """An output from an older version has unknown parameters, not wrong ones."""
    output_path, organism_name = organism
    write_roary_inputs(output_path, organism_name, {GENES[0]: 10}, genomes=10)

    genes_csv = os.path.join(
        output_path, organism_name, 'conservation', 'roary_output', 'results',
        'gene_presence_absence.csv'
    )
    files.write_manifest(genes_csv, {'min_core_freq': 99, 'min_identity': 95})

    genome.core_genome_programs(output_path, organism_name, core_threshold=99, identity=95,
                                cpus=1, program_list=['roary'])
    captured = capsys.readouterr().out
    assert 'genome_count' in captured and 'unknown' in captured
    assert os.path.exists(genes_csv)
