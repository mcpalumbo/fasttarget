"""Microbiome plot helpers: file names, percentages and category order.

These build the per-gene plot files and the composition tables, so a name that
escapes its directory or a denominator of zero would show up as a crash in the
middle of a long run, or as a file written somewhere unintended.

Ported from an earlier local suite and extended; this module had no coverage in
the fast suite before.
"""

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ftscripts import microbiome_hit_plots as plots


# -------------------------------------------------------------- file names

@pytest.mark.parametrize('value,expected', [
    (' ../gene A/B ', 'gene_A_B'),
    ('MPN_RS00140', 'MPN_RS00140'),
    ('a b  c', 'a_b_c'),
    ('a\tb', 'a_b'),
    # Dots survive, separators do not: the result is still a single component.
    ('gene/../../etc/passwd', 'gene_.._.._etc_passwd'),
])
def test_a_name_is_reduced_to_a_single_safe_component(value, expected):
    """
    The locus tag becomes a file name, so a separator or a parent reference in it
    would write outside the plot directory.
    """
    safe = plots.safe_filename(value)

    assert safe == expected
    assert '/' not in safe and '\\' not in safe
    assert safe not in {'.', '..'}


@pytest.mark.parametrize('value', ['...', '.', '..', '', '   ', '///'])
def test_a_name_with_nothing_usable_falls_back_to_a_fixed_label(value):
    assert plots.safe_filename(value) == 'locus_tag'


# ------------------------------------------------------------- percentages

@pytest.mark.parametrize('part,whole,expected', [
    (1, 4, 25.0),
    (0, 10, 0.0),
    (10, 10, 100.0),
    (1, 3, pytest.approx(33.333, abs=0.001)),
])
def test_a_percentage_is_the_part_over_the_whole(part, whole, expected):
    assert plots.percentage(part, whole) == expected


@pytest.mark.parametrize('whole', [0, None])
def test_an_empty_denominator_is_zero_rather_than_an_error(whole):
    """
    A rank with no classified groups has no denominator. Zero is the readable
    answer; dividing would end the plotting stage mid-run.
    """
    assert plots.percentage(1, whole) == 0.0


# --------------------------------------------------------- category order

def test_other_and_unclassified_are_sorted_last():
    """
    They are not taxa, they are the remainder, so they belong at the end of a
    composition bar however the real names sort.
    """
    ordered = plots.order_composition_columns(
        [plots.OTHER, 'Beta', plots.UNCLASSIFIED, 'Alpha']
    )

    assert ordered == ['Alpha', 'Beta', plots.UNCLASSIFIED, plots.OTHER]


def test_real_categories_keep_alphabetical_order():
    ordered = plots.order_composition_columns(['Zeta', 'Alpha', 'Mu'])

    assert ordered == ['Alpha', 'Mu', 'Zeta']


def test_the_remainder_categories_are_only_added_when_present():
    assert plots.order_composition_columns(['Alpha']) == ['Alpha']
    assert plots.order_composition_columns([plots.OTHER]) == [plots.OTHER]
    assert plots.order_composition_columns([]) == []


# ------------------------------------------------------- parquet validation

def test_a_parquet_with_the_expected_columns_is_accepted(tmp_path):
    path = tmp_path / 'data.parquet'
    pq.write_table(pa.table({'gene': ['g1']}), path)

    plots.validate_parquet_columns(path, ['gene'], 'Data')


def test_a_parquet_missing_a_column_is_reported_by_name(tmp_path):
    """
    The plots read columns the consolidation step produced; a schema change has
    to say which column is gone rather than fail later on a KeyError.
    """
    path = tmp_path / 'data.parquet'
    pq.write_table(pa.table({'gene': ['g1']}), path)

    with pytest.raises(SystemExit, match='missing expected columns'):
        plots.validate_parquet_columns(path, ['gene', 'pident'], 'Data')
