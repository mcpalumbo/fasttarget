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


def test_none_tables_are_ignored(organism):
    output_path, organism_name = organism
    fasttarget.reset_failed_steps()

    good = _module_table(GENES, 'ok', [True] * len(GENES))
    results = fasttarget.merge_final_tables(StubConfig(organism_name), output_path, [None, good, None])

    assert len(results) == len(GENES)
    assert 'ok' in results.columns
