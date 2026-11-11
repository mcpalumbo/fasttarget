"""Sequence identity, BLAST cutoffs and the UniParc tier ordering.

These are the decisions that determine whether a structure is assigned to a
gene, and the requirement they implement is that the result describes the user's
own protein: an exact sequence match, not a homolog and not a domain.
"""

import os

import pytest

from conftest import GENES, SEQUENCES, write_blast_output
from ftscripts import structures


# ------------------------------------------------- sequence normalisation

def test_whitespace_and_a_trailing_stop_are_removed():
    assert structures.normalize_protein_sequence('MKA ILV\nVLL\t') == 'MKAILVVLL'
    assert structures.normalize_protein_sequence('MKAILVVLL*') == 'MKAILVVLL'


def test_case_is_normalised_but_residues_are_never_dropped():
    """
    Only whitespace and the stop codon go: the route depends on byte-level
    identity, so an unusual residue must not be silently cleaned away.
    """
    assert structures.normalize_protein_sequence('mkailvvll') == 'MKAILVVLL'

    for residue in 'XUBZJO':
        sequence = f'MKA{residue}ILV'
        assert structures.normalize_protein_sequence(sequence) == sequence


def test_an_internal_stop_is_kept():
    """A trailing stop is formatting; an internal one changes the protein."""
    assert structures.normalize_protein_sequence('MKA*ILV') == 'MKA*ILV'
    assert structures.normalize_protein_sequence('MKA*ILV*') == 'MKA*ILV'


# ---------------------------------------------------------------- checksum

def _reflected_crc64(sequence, polynomial=0xd800000000000000):
    """
    The SWISS-PROT reflected CRC64, written out from its polynomial.

    Independent of the library the pipeline calls, so this pins the variant: the
    checksum is the key UniParc is indexed on, and computing the wrong one would
    not fail loudly, it would simply never find a match.
    """
    crc = 0
    for byte in sequence.encode('ascii'):
        crc ^= byte
        for _ in range(8):
            crc = (crc >> 1) ^ (polynomial if crc & 1 else 0)
    return f'{crc:016X}'


def _iso_crc64(sequence, polynomial=0x1B):
    """CRC64-ISO, the variant this is not."""
    crc = 0
    mask = 0xFFFFFFFFFFFFFFFF
    for byte in sequence.encode('ascii'):
        crc ^= byte << 56
        for _ in range(8):
            crc = ((crc << 1) ^ polynomial) & mask if crc & (1 << 63) else (crc << 1) & mask
    return f'{crc:016X}'


def test_the_checksum_is_the_variant_uniparc_indexes_on():
    insulin_a_chain = 'GIVEQCCTSICSLYQLENYCN'
    checksum = structures.protein_sequence_checksum(insulin_a_chain)

    assert checksum == 'AC2453F8AFA43016'
    assert checksum == _reflected_crc64(insulin_a_chain)
    assert checksum != _iso_crc64(insulin_a_chain), 'CRC64-ISO would never match UniParc'
    assert len(checksum) == 16 and 'CRC-' not in checksum


def test_the_checksum_ignores_formatting_like_uniparc_does():
    insulin_a_chain = 'GIVEQCCTSICSLYQLENYCN'
    assert structures.protein_sequence_checksum(insulin_a_chain) == \
        structures.protein_sequence_checksum(f'  {insulin_a_chain.lower()}  *')


def test_different_sequences_get_different_checksums():
    checksums = {structures.protein_sequence_checksum(seq) for seq in SEQUENCES.values()}
    assert len(checksums) == len(SEQUENCES)


def test_a_single_residue_change_changes_the_checksum():
    """What makes the tier safe: the checksum cannot collapse two proteins."""
    first = structures.protein_sequence_checksum('MKAILVVLLYTFATANA')
    second = structures.protein_sequence_checksum('MKAILVVLLYTFATANG')
    assert first != second


# ------------------------------------------------------------ BLAST cutoffs

def _blast_file(tmp_path, rows):
    """Write a BLAST file with qlen and slen, which this parser reads."""
    path = os.path.join(str(tmp_path), 'blast.tsv')
    with open(path, 'w') as handle:
        for qseqid, sseqid, pident, qcovs, length, slen in rows:
            handle.write('\t'.join(str(v) for v in [
                qseqid, sseqid, pident, length, 0, 0, 1, length, 1, length,
                '1e-50', 200.0, qcovs, qcovs, length, slen,
            ]) + '\n')
    return path


def test_a_hit_passing_every_cutoff_is_returned(tmp_path):
    path = _blast_file(tmp_path, [(GENES[0], 'P12345', 99.0, 100, 100, 100)])
    assert structures.parse_best_result_blast(path) == {GENES[0]: 'P12345'}


def test_cutoffs_are_inclusive_so_a_value_exactly_at_the_limit_passes(tmp_path):
    """The comparison is '>=', which only a boundary case distinguishes."""
    path = _blast_file(tmp_path, [(GENES[0], 'P12345', 95.0, 90, 90, 100)])
    assert structures.parse_best_result_blast(path) == {GENES[0]: 'P12345'}


@pytest.mark.parametrize('pident,qcovs,length,slen,reason', [
    (94.9, 100, 100, 100, 'identity just below the cutoff'),
    (99.0, 89, 89, 100, 'query coverage just below the cutoff'),
    (99.0, 100, 50, 100, 'subject coverage below the cutoff: a domain-only match'),
])
def test_a_hit_failing_any_cutoff_is_not_a_mapping(tmp_path, pident, qcovs, length, slen, reason):
    path = _blast_file(tmp_path, [(GENES[0], 'P12345', pident, qcovs, length, slen)])
    assert structures.parse_best_result_blast(path) == {}, reason


def test_a_domain_match_over_a_long_subject_is_rejected(tmp_path):
    """
    Full query coverage against a much longer subject means the gene matches one
    domain of a larger protein, which is not the user's protein.
    """
    path = _blast_file(tmp_path, [(GENES[0], 'P12345', 100.0, 100, 100, 1000)])
    assert structures.parse_best_result_blast(path) == {}


def test_the_best_hit_is_the_highest_identity_then_coverage(tmp_path):
    path = _blast_file(tmp_path, [
        (GENES[0], 'LOW', 96.0, 100, 100, 100),
        (GENES[0], 'BEST', 99.0, 100, 100, 100),
        (GENES[0], 'MID', 97.0, 100, 100, 100),
    ])
    assert structures.parse_best_result_blast(path) == {GENES[0]: 'BEST'}


def test_all_hits_returns_every_passing_subject(tmp_path):
    path = _blast_file(tmp_path, [
        (GENES[0], 'A', 99.0, 100, 100, 100),
        (GENES[0], 'B', 99.0, 100, 100, 100),
        (GENES[0], 'C', 96.0, 100, 100, 100),
    ])
    result = structures.parse_best_result_blast(path, all_hits=True)
    assert sorted(result[GENES[0]]) == ['A', 'B']


def test_an_empty_or_absent_blast_file_is_an_empty_result(tmp_path):
    missing = os.path.join(str(tmp_path), 'absent.tsv')
    assert structures.parse_best_result_blast(missing) == {}

    empty = os.path.join(str(tmp_path), 'empty.tsv')
    open(empty, 'w').close()
    assert structures.parse_best_result_blast(empty) == {}


def test_each_query_is_mapped_independently(tmp_path):
    path = _blast_file(tmp_path, [
        (GENES[0], 'P1', 99.0, 100, 100, 100),
        (GENES[1], 'P2', 50.0, 100, 100, 100),
        (GENES[2], 'P3', 99.0, 100, 100, 100),
    ])
    result = structures.parse_best_result_blast(path)
    assert result == {GENES[0]: 'P1', GENES[2]: 'P3'}


# ------------------------------------------------------- UniParc tier order

def _xref(accession, taxid, active=True, database='UniProtKB/TrEMBL'):
    return {
        'database': database,
        'id': accession,
        'active': active,
        'organism': {'taxonId': taxid},
    }


SPECIES, STRAIN, OTHER = 2104, 272634, 9606


def test_the_strain_outranks_the_species_and_anything_else():
    entry = {'uniParcCrossReferences': [
        _xref('OTHER_ORG', OTHER),
        _xref('SPECIES', SPECIES),
        _xref('STRAIN', STRAIN),
    ]}

    ranked = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)

    assert [c['accession'] for c in ranked] == ['STRAIN', 'SPECIES', 'OTHER_ORG']
    assert [c['tier'] for c in ranked] == [4, 5, 6]


def test_an_active_entry_outranks_a_withdrawn_one_in_the_same_tier():
    entry = {'uniParcCrossReferences': [
        _xref('WITHDRAWN', STRAIN, active=False),
        _xref('ACTIVE', STRAIN, active=True),
    ]}

    ranked = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)

    assert [c['accession'] for c in ranked] == ['ACTIVE', 'WITHDRAWN']


def test_a_withdrawn_strain_entry_still_outranks_an_active_foreign_one():
    """Taxonomic closeness is the first key; being active only breaks ties."""
    entry = {'uniParcCrossReferences': [
        _xref('FOREIGN_ACTIVE', OTHER, active=True),
        _xref('STRAIN_WITHDRAWN', STRAIN, active=False),
    ]}

    ranked = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)

    assert [c['accession'] for c in ranked] == ['STRAIN_WITHDRAWN', 'FOREIGN_ACTIVE']


def test_ordering_is_deterministic_for_equal_candidates():
    entry = {'uniParcCrossReferences': [
        _xref('B00002', STRAIN), _xref('A00001', STRAIN), _xref('C00003', STRAIN),
    ]}

    ranked = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)

    assert [c['accession'] for c in ranked] == ['A00001', 'B00002', 'C00003']


def test_non_uniprot_cross_references_are_ignored():
    """A UniParc entry also lists EMBL, RefSeq and patent records."""
    entry = {'uniParcCrossReferences': [
        _xref('EMBL_ID', STRAIN, database='EMBL'),
        _xref('REFSEQ_ID', STRAIN, database='RefSeq'),
        {'database': 'EMBLWGS', 'id': 'WGS_ID', 'active': True},
        _xref('KB_ID', STRAIN, database='UniProtKB/Swiss-Prot'),
    ]}

    ranked = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)

    assert [c['accession'] for c in ranked] == ['KB_ID']


def test_a_cross_reference_with_no_accession_is_skipped():
    entry = {'uniParcCrossReferences': [
        {'database': 'UniProtKB/TrEMBL', 'active': True, 'organism': {'taxonId': STRAIN}},
        _xref('GOOD', STRAIN),
    ]}

    ranked = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)

    assert [c['accession'] for c in ranked] == ['GOOD']


def test_an_entry_without_cross_references_ranks_nothing():
    assert structures.rank_uniparc_uniprot_xrefs({}, SPECIES, STRAIN) == []
    assert structures.rank_uniparc_uniprot_xrefs({'uniParcCrossReferences': None}, SPECIES, STRAIN) == []


def test_taxids_compare_as_text_so_numbers_and_strings_agree():
    """The API returns integers; the config may supply either."""
    entry = {'uniParcCrossReferences': [_xref('STRAIN', STRAIN)]}

    as_numbers = structures.rank_uniparc_uniprot_xrefs(entry, SPECIES, STRAIN)
    as_text = structures.rank_uniparc_uniprot_xrefs(entry, str(SPECIES), str(STRAIN))

    assert as_numbers[0]['tier'] == as_text[0]['tier'] == 4
