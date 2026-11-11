"""BindingDB: dump parsing, measurement interpretation and hit selection.

Two hand-written parsers carry the risk here. The MySQL literal reader turns the
official dump into rows, so a mishandled escape corrupts data silently. And
parse_measurement decides whether an affinity is an exact number or only a
bound, which is what makes a compound's potency comparable or not.
"""

import math

import pandas as pd
import pytest

from ftscripts import bindingdb


# ------------------------------------------------------ MySQL dump literals

def _rows(values):
    return list(bindingdb.iter_mysql_rows(values))


def test_a_simple_row_is_decoded():
    assert _rows("(1,'abc',2.5)") == [['1', 'abc', '2.5']]


def test_several_rows_are_decoded_in_order():
    assert _rows("(1,'a'),(2,'b'),(3,'c')") == [['1', 'a'], ['2', 'b'], ['3', 'c']]


def test_a_trailing_semicolon_ends_the_statement():
    assert _rows("(1,'a'),(2,'b');") == [['1', 'a'], ['2', 'b']]


def test_null_becomes_none_and_is_not_the_text_null():
    row = _rows("(NULL,'NULL')")[0]

    assert row[0] is None
    assert row[1] == 'NULL', 'a quoted NULL is the string'


def test_a_quoted_number_stays_a_string():
    """
    Identifiers such as a polymer id arrive quoted; converting them to numbers
    would lose leading zeros and change how they join.
    """
    row = _rows("('007',7)")[0]

    assert row[0] == '007'
    assert row[1] == '7'


@pytest.mark.parametrize('literal,expected', [
    (r"('a\nb')", 'a\nb'),
    (r"('a\tb')", 'a\tb'),
    (r"('a\rb')", 'a\rb'),
    (r"('a\0b')", 'a\0b'),
    (r"('a\\b')", 'a\\b'),
    (r"('a\'b')", "a'b"),
    (r"('a\Zb')", 'a\x1ab'),
    (r"('a\qb')", 'aqb'),
])
def test_backslash_escapes_are_decoded(literal, expected):
    assert _rows(literal)[0][0] == expected


def test_a_doubled_quote_is_one_quote():
    assert _rows("('it''s')")[0][0] == "it's"


def test_unicode_survives_decoding():
    """Compound and organism names carry accents and Greek letters."""
    assert _rows("('β-lactam, Müller')")[0][0] == 'β-lactam, Müller'


def test_a_comma_inside_a_quoted_value_is_not_a_separator():
    assert _rows("('a,b','c')")[0] == ['a,b', 'c']


def test_a_parenthesis_inside_a_quoted_value_does_not_end_the_row():
    assert _rows("('f(x)','y')")[0] == ['f(x)', 'y']


def test_whitespace_between_tokens_is_ignored():
    assert _rows(" ( 1 , 'a' ) , ( 2 , 'b' ) ") == [['1', 'a'], ['2', 'b']]


def test_scientific_notation_and_signs_are_kept_verbatim():
    row = _rows("(1e-9,-2.5,+3,.5)")[0]

    assert row == ['1e-9', '-2.5', '+3', '.5']


@pytest.mark.parametrize('bad,fragment', [
    ("1,'a'", 'Expected a row'),
    ("(1,'a'", 'Unterminated'),
    ("(1 'a')", 'value separator'),
    ("(1,'a'),", 'Trailing comma'),
    ("(1,@)", 'Unsupported MySQL literal'),
    ("(1,'a') garbage", 'Unexpected SQL'),
])
def test_malformed_input_is_refused_rather_than_guessed(bad, fragment):
    """
    A dump this parser cannot read must stop the preparation: a row decoded
    wrongly would be indistinguishable from real data afterwards.
    """
    with pytest.raises(ValueError, match=fragment):
        _rows(bad)


def test_an_empty_statement_yields_nothing():
    assert _rows('') == []
    assert _rows('   ') == []
    assert _rows(';') == []


# ------------------------------------------------------------- measurements

def test_an_exact_value_gets_a_pactivity():
    """pActivity is only meaningful for an exact concentration."""
    result = bindingdb.parse_measurement('100')

    assert result['value'] == 100.0
    assert result['relation'] == '='
    assert result['value_status'] == 'exact'
    assert result['pactivity'] == pytest.approx(9 - math.log10(100))
    assert result['units'] == 'nM'


@pytest.mark.parametrize('raw,relation', [
    ('>100', '>'), ('<100', '<'), ('>=100', '>='), ('<=100', '<='), ('~100', '~'),
])
def test_a_censored_value_is_kept_as_a_bound_without_pactivity(raw, relation):
    """
    A measurement reported as '>10000 nM' is an upper limit on potency, not a
    potency: giving it a pActivity would make it comparable to exact values.
    """
    result = bindingdb.parse_measurement(raw)

    assert result['relation'] == relation
    assert result['value'] == 100.0
    assert result['value_status'] == 'bound_or_approximate'
    assert result['pactivity'] is None


@pytest.mark.parametrize('raw,relation', [('≤100', '<='), ('≥100', '>=')])
def test_unicode_comparison_signs_are_understood(raw, relation):
    assert bindingdb.parse_measurement(raw)['relation'] == relation


def test_html_entities_are_decoded_before_parsing():
    """The dump carries HTML-escaped text in places."""
    assert bindingdb.parse_measurement('&gt;100')['relation'] == '>'


@pytest.mark.parametrize('raw', [None, '', '   '])
def test_a_missing_value_is_recorded_as_missing(raw):
    result = bindingdb.parse_measurement(raw)

    assert result['value_status'] == 'missing'
    assert result['value'] is None


@pytest.mark.parametrize('raw', ['n/a', 'abc', '10 to 20', '1,000'])
def test_an_unreadable_value_stays_unparsed_rather_than_becoming_a_number(raw):
    result = bindingdb.parse_measurement(raw)

    assert result['value_status'] == 'unparsed'
    assert result['value'] is None
    assert result['value_raw'] == raw, 'the original is kept for inspection'


@pytest.mark.parametrize('raw', ['0', '-5', '0.0'])
def test_a_non_positive_concentration_is_rejected(raw):
    """
    A concentration of zero or less is not measurable, and log10 of it is not a
    number: without this check pActivity would be infinite or complex.
    """
    result = bindingdb.parse_measurement(raw)

    assert result['value_status'] == 'nonpositive_or_nonfinite'
    assert result['pactivity'] is None


def test_a_relation_that_contradicts_the_qualifier_is_refused():
    """The row says both '>' and '<', so neither can be trusted."""
    result = bindingdb.parse_measurement('>100', uncertainty='<')

    assert result['value_status'] == 'conflicting_relations'
    assert result['value'] is None


def test_a_matching_qualifier_is_accepted():
    result = bindingdb.parse_measurement('>100', uncertainty='>')

    assert result['relation'] == '>'
    assert result['value'] == 100.0


def test_a_qualifier_alone_supplies_the_relation():
    result = bindingdb.parse_measurement('100', uncertainty='>')

    assert result['relation'] == '>'
    assert result['pactivity'] is None


def test_an_unsupported_qualifier_keeps_the_value_but_flags_it():
    result = bindingdb.parse_measurement('100', uncertainty='approximately')

    assert result['value'] == 100.0
    assert result['value_status'] == 'unsupported_uncertainty'
    assert result['pactivity'] is None


def test_pactivity_follows_the_standard_conversion():
    """1 nM is pActivity 9, 1000 nM is 6."""
    assert bindingdb.parse_measurement('1')['pactivity'] == pytest.approx(9.0)
    assert bindingdb.parse_measurement('1000')['pactivity'] == pytest.approx(6.0)


# ------------------------------------------------------- protein evidence

@pytest.mark.parametrize('value,expected', [
    ('MKAILV', 'MKAILV'),
    ('mka ilv*', 'MKAILV'),
    (None, ''),
    ('', ''),
])
def test_sequences_are_normalised(value, expected):
    assert bindingdb._normalized_sequence(value) == expected


def test_an_amino_acid_sequence_is_usable_evidence():
    assert bindingdb._protein_evidence({'sequence': 'MKAILVVLLYTFATANA'}) == 'amino_acid_sequence'


@pytest.mark.parametrize('sequence', ['ACGTACGTACGT', 'ACGUACGU', 'NNNNNNN', 'ACGTN'])
def test_a_nucleotide_sequence_is_not_protein_evidence(sequence):
    """
    A polymer whose sequence uses only nucleotide letters is DNA or RNA; aligning
    the genome's proteins against it would be meaningless.
    """
    assert bindingdb._protein_evidence({'sequence': sequence}) is None


@pytest.mark.parametrize('sequence', ['', None, 'MKA-ILV', 'MKA ILV 123', 'MKA@ILV'])
def test_an_unusable_sequence_is_not_evidence(sequence):
    assert bindingdb._protein_evidence({'sequence': sequence}) is None


def test_the_free_text_type_does_not_override_the_sequence():
    """Only the sequence decides, so a mislabelled row cannot slip through."""
    row = {'sequence': 'ACGTACGT', 'target_type': 'protein', 'unpid1': 'P12345'}

    assert bindingdb._protein_evidence(row) is None


# ---------------------------------------------------------- best hit choice

def _candidates(rows):
    return pd.DataFrame(rows, columns=['polymerid', 'gene', 'bitscore', 'evalue'])


def test_the_highest_bitscore_wins():
    candidates = _candidates([
        (1, 'GENE_A', 100.0, 1e-30),
        (1, 'GENE_B', 200.0, 1e-50),
    ])

    best, ambiguous = bindingdb.select_best_hits(candidates, 'unassigned')

    assert list(best.gene) == ['GENE_B']
    assert ambiguous.empty


def test_the_lowest_evalue_breaks_a_bitscore_tie():
    candidates = _candidates([
        (1, 'GENE_A', 200.0, 1e-30),
        (1, 'GENE_B', 200.0, 1e-50),
    ])

    best, ambiguous = bindingdb.select_best_hits(candidates, 'unassigned')

    assert list(best.gene) == ['GENE_B']
    assert ambiguous.empty


def test_an_exact_tie_is_left_unassigned_by_default():
    """
    Two genes equally good for one BindingDB polymer is a real ambiguity;
    picking one would attribute a compound to an arbitrary gene.
    """
    candidates = _candidates([
        (1, 'GENE_A', 200.0, 1e-50),
        (1, 'GENE_B', 200.0, 1e-50),
    ])

    best, ambiguous = bindingdb.select_best_hits(candidates, 'unassigned')

    assert best.empty
    assert sorted(ambiguous.gene) == ['GENE_A', 'GENE_B']


def test_a_tie_can_be_kept_for_every_gene_on_request():
    candidates = _candidates([
        (1, 'GENE_A', 200.0, 1e-50),
        (1, 'GENE_B', 200.0, 1e-50),
    ])

    best, ambiguous = bindingdb.select_best_hits(candidates, 'all')

    assert sorted(best.gene) == ['GENE_A', 'GENE_B']
    assert sorted(ambiguous.gene) == ['GENE_A', 'GENE_B']


def test_each_polymer_is_resolved_independently():
    candidates = _candidates([
        (1, 'GENE_A', 200.0, 1e-50),
        (2, 'GENE_B', 100.0, 1e-20),
        (2, 'GENE_C', 100.0, 1e-20),
    ])

    best, ambiguous = bindingdb.select_best_hits(candidates, 'unassigned')

    assert list(best.polymerid) == [1]
    assert sorted(ambiguous.polymerid) == [2, 2]


def test_a_repeated_pair_is_counted_once():
    candidates = _candidates([
        (1, 'GENE_A', 200.0, 1e-50),
        (1, 'GENE_A', 200.0, 1e-50),
    ])

    best, ambiguous = bindingdb.select_best_hits(candidates, 'unassigned')

    assert list(best.gene) == ['GENE_A'], 'a duplicate row is not a tie'
    assert ambiguous.empty


def test_an_unknown_tie_policy_is_refused():
    with pytest.raises(ValueError, match='tie_policy'):
        bindingdb.select_best_hits(_candidates([]), 'pick_one')


def test_no_candidates_yields_nothing():
    best, ambiguous = bindingdb.select_best_hits(_candidates([]), 'unassigned')

    assert best.empty and ambiguous.empty


# ------------------------------------------------------------ mapping guards

def _fasta(tmp_path, records):
    path = tmp_path / 'proteins.faa'
    path.write_text(''.join(f'>{name}\n{sequence}\n' for name, sequence in records))
    return path


@pytest.mark.parametrize('identity,coverage,threads', [
    (0, 70, 4), (101, 70, 4), (90, 0, 4), (90, 101, 4), (90, 70, 0),
])
def test_out_of_range_mapping_parameters_are_refused(tmp_path, identity, coverage, threads):
    fasta = _fasta(tmp_path, [('GENE1', 'MKAILV')])
    targets = pd.DataFrame(columns=bindingdb.PROTEIN_COLUMNS)

    with pytest.raises(ValueError, match='Identity/coverage'):
        bindingdb.map_proteins(fasta, targets, tmp_path / 'work',
                               min_identity=identity, min_coverage=coverage, threads=threads)
