"""The merged results table.

Every join is a left join, so whichever table comes first decides the gene set.
Anchoring on the annotated CDS list is what makes the output one row per gene
regardless of which modules ran, in what order.
"""

import os

import pandas as pd
import pytest

import configuration
import fasttarget
from conftest import GENES


class StubConfig:
    """The attributes merge_final_tables and the summary read."""

    def __init__(self, organism_name):
        self.organism = {'name': organism_name, 'tax_id': 1, 'strain_taxid': 2}
        self.metabolism_pathwaytools = None
        self.metabolism_sbml = None
        self.structures = None
        self.ligands = None
        self.core = None
        self.offtarget = None
        self.deg = None
        self.psortb = None
        self.metadata = None
        self.colabfold = None


def _module_table(genes, column, values):
    return pd.DataFrame({'gene': genes, column: values})


def test_every_annotated_gene_has_a_row_whatever_the_modules_returned(organism):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()

    # A module that only knows about two genes must not shrink the table.
    partial = _module_table(GENES[:2], 'core_roary', [True, False])

    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [partial])

    assert sorted(results['gene']) == sorted(GENES)
    assert len(results) == len(GENES)


def test_module_order_does_not_change_the_gene_set(organism):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()

    first = _module_table(GENES[:2], 'a', [1, 2])
    second = _module_table(GENES, 'b', list(range(len(GENES))))

    one_way = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [first, second])
    other_way = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [second, first])

    assert sorted(one_way['gene']) == sorted(other_way['gene']) == sorted(GENES)


def test_a_gene_absent_from_the_annotation_is_dropped_with_a_warning(organism, caplog):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()

    stray = _module_table(GENES + ['NOT_IN_ANNOTATION'], 'x', [1] * (len(GENES) + 1))

    with caplog.at_level('WARNING'):
        results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [stray])

    assert 'NOT_IN_ANNOTATION' not in set(results['gene'])
    assert any('absent from the reference annotation' in record.message for record in caplog.records)


def test_no_module_tables_produces_no_results_rather_than_a_gene_list(organism):
    """
    The anchor is not a result: a run where every module failed must not look
    like a successful run that produced a list of genes.
    """
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()

    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [])

    assert results is None


def test_a_repeated_gene_does_not_multiply_rows(organism):
    """
    Every join is a left join on 'gene', so a module table carrying a gene twice
    multiplied that gene's rows and the table stopped having one row per gene.
    """
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()
    repeated = pd.DataFrame({
        'gene': [GENES[0], GENES[0], GENES[1]],
        'core_roary': [True, False, True],
    })

    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [repeated])

    assert len(results) == len(GENES)
    assert (results['gene'] == GENES[0]).sum() == 1
    assert results.loc[results['gene'] == GENES[0], 'core_roary'].iloc[0] is True, 'the first row'


def test_a_repeated_gene_is_reported_as_a_failure(organism):
    """
    The two rows disagreed, so one of them is wrong: dropping it quietly would
    hide a module defect.
    """
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()
    repeated = pd.DataFrame({'gene': [GENES[0], GENES[0]], 'core_roary': [True, False]})

    fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [repeated])

    steps = fasttarget.failed_steps()
    assert steps, 'a repeated gene must make the run exit non-zero'
    assert 'core_roary' in steps[0], 'the offending table has to be identifiable'


def test_a_column_supplied_twice_keeps_its_documented_name(organism):
    """
    pandas renames both copies to score_x and score_y without an error, so
    neither matches the name the module documented.
    """
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()
    first = _module_table(GENES, 'score', list(range(len(GENES))))
    second = _module_table(GENES, 'score', [99] * len(GENES))

    results = fasttarget.merge_final_tables(
        StubConfig(organism_name), output_path, [first, second]
    )

    assert 'score' in results.columns
    assert not [name for name in results.columns if name.endswith(('_x', '_y'))]
    assert results['score'].tolist() == list(range(len(GENES))), 'the first module wins'
    assert any('score' in step for step in fasttarget.failed_steps())


def test_a_clean_merge_reports_nothing(organism):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()
    clean = _module_table(GENES, 'core_roary', [True] * len(GENES))

    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [clean])

    assert len(results) == len(GENES)
    assert fasttarget.failed_steps() == []


def test_distinct_columns_from_several_modules_are_all_kept(organism):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()
    tables = [
        _module_table(GENES, 'core_roary', [True] * len(GENES)),
        _module_table(GENES, 'hit_in_deg', [False] * len(GENES)),
        _module_table(GENES, 'druggability_score', [0.5] * len(GENES)),
    ]

    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, tables)

    for column in ('core_roary', 'hit_in_deg', 'druggability_score'):
        assert column in results.columns
    assert fasttarget.failed_steps() == []


def test_none_tables_are_ignored(organism):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()

    good = _module_table(GENES, 'ok', [True] * len(GENES))
    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [None, good, None])

    assert len(results) == len(GENES)
    assert 'ok' in results.columns
