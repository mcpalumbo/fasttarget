"""ChEMBL: schema guard, activity filters and evidence labelling.

The filters decide which measurements count as evidence that a compound acts on
a target, so each reason has to apply on its own and the audit has to keep what
was excluded and why. Nothing here opens the real ChEMBL database: the schema
tests build the required tables in memory, and the rest work on DataFrames.
"""

import sqlite3

import pandas as pd
import pytest

from ftscripts import chembl


# -------------------------------------------------------------- schema guard

def _chembl_schema(connection, release='ChEMBL_35', omit=None):
    """Create the tables and columns the module declares as required."""
    for table, columns in chembl.REQUIRED.items():
        wanted = set(columns)
        if omit and omit[0] == table:
            wanted = wanted - {omit[1]}
        definition = ', '.join(f'"{column}"' for column in sorted(wanted))
        connection.execute(f'CREATE TABLE "{table}" ({definition})')
    if release:
        connection.execute('INSERT INTO version (name) VALUES (?)', (release,))
    connection.commit()


@pytest.fixture
def chembl_db():
    connection = sqlite3.connect(':memory:')
    yield connection
    connection.close()


def test_a_supported_schema_is_accepted_with_its_release(chembl_db):
    _chembl_schema(chembl_db)

    schema, release = chembl.inspect_schema(chembl_db)

    assert set(schema) == set(chembl.REQUIRED)
    assert release == ['ChEMBL_35']


@pytest.mark.parametrize('table,column', [
    ('target_dictionary', 'tax_id'),
    ('activities', 'pchembl_value'),
    ('assays', 'confidence_score'),
    ('component_sequences', 'sequence'),
])
def test_a_missing_column_is_refused_naming_table_and_column(chembl_db, table, column):
    """
    A new ChEMBL release that drops or renames a column must stop the stage
    rather than produce a silently different answer.
    """
    _chembl_schema(chembl_db, omit=(table, column))

    with pytest.raises(ValueError) as error:
        chembl.inspect_schema(chembl_db)

    message = str(error.value)
    assert table in message and column in message


def test_a_database_without_a_release_identifier_is_refused(chembl_db):
    _chembl_schema(chembl_db, release=None)

    with pytest.raises(ValueError, match='release identifier'):
        chembl.inspect_schema(chembl_db)


def test_a_version_row_that_is_not_a_chembl_release_does_not_count(chembl_db):
    _chembl_schema(chembl_db, release='something else')

    with pytest.raises(ValueError, match='release identifier'):
        chembl.inspect_schema(chembl_db)


# ------------------------------------------------------ sequence handling

@pytest.mark.parametrize('value,expected', [
    ('MKAILV', 'MKAILV'),
    ('mkailv', 'MKAILV'),
    ('MKA ILV\n', 'MKAILV'),
    ('MKAILV*', 'MKAILV'),
    ('  mka\tilv*  ', 'MKAILV'),
    (None, ''),
    (float('nan'), ''),
])
def test_sequences_are_normalised_for_comparison(value, expected):
    assert chembl._sequence(value) == expected


# ------------------------------------------------------- activity filtering

def _activities(**overrides):
    row = {
        'activity_id': 1,
        'relationship_type': 'D',
        'confidence_score': 9,
        'pchembl_value': 7.5,
        'data_validity_comment': None,
    }
    row.update(overrides)
    return pd.DataFrame([row])


def test_a_clean_measurement_is_retained():
    audit, retained = chembl.filter_activities(_activities())

    assert list(audit.filter_status) == ['retained']
    assert len(retained) == 1
    assert list(audit.filter_reasons) == ['']


def test_an_indirect_assignment_is_excluded_only_when_asked():
    frame = _activities(relationship_type='H')

    _, kept = chembl.filter_activities(frame, direct_only=False)
    assert len(kept) == 1

    audit, dropped = chembl.filter_activities(frame, direct_only=True)
    assert len(dropped) == 0
    assert 'not_direct_assignment' in audit.filter_reasons.iloc[0]


@pytest.mark.parametrize('score,minimum,excluded', [
    (9, 8, False),
    (8, 8, False),      # the minimum is inclusive
    (7, 8, True),
    (None, 8, True),    # an unknown score cannot satisfy a minimum
])
def test_the_confidence_minimum_is_inclusive(score, minimum, excluded):
    frame = _activities(confidence_score=score)

    audit, retained = chembl.filter_activities(frame, min_confidence_score=minimum)

    assert (len(retained) == 0) is excluded
    if excluded:
        assert 'confidence_below_minimum' in audit.filter_reasons.iloc[0]


def test_a_missing_pchembl_is_excluded_only_when_required():
    frame = _activities(pchembl_value=None)

    _, kept = chembl.filter_activities(frame, require_pchembl=False)
    assert len(kept) == 1

    audit, dropped = chembl.filter_activities(frame, require_pchembl=True)
    assert len(dropped) == 0
    assert 'missing_pchembl' in audit.filter_reasons.iloc[0]


@pytest.mark.parametrize('comment,excluded', [
    (None, False),
    ('Manually validated', False),   # a comment that is not a problem
    ('Outside typical range', True),
    ('Potential author error', True),
])
def test_a_data_validity_comment_excludes_unless_manually_validated(comment, excluded):
    frame = _activities(data_validity_comment=comment)

    audit, retained = chembl.filter_activities(frame, exclude_problematic=True)

    assert (len(retained) == 0) is excluded
    if excluded:
        assert 'problematic_data_validity' in audit.filter_reasons.iloc[0]


def test_problematic_data_can_be_kept_deliberately():
    frame = _activities(data_validity_comment='Outside typical range')

    _, retained = chembl.filter_activities(frame, exclude_problematic=False)

    assert len(retained) == 1


def test_every_reason_for_exclusion_is_recorded_not_just_the_first():
    """The audit explains the decision, so a row can carry several reasons."""
    frame = _activities(relationship_type='H', confidence_score=3, pchembl_value=None,
                        data_validity_comment='Outside typical range')

    audit, retained = chembl.filter_activities(
        frame, direct_only=True, min_confidence_score=8,
        require_pchembl=True, exclude_problematic=True,
    )

    reasons = audit.filter_reasons.iloc[0].split(';')
    assert set(reasons) == {
        'not_direct_assignment', 'confidence_below_minimum',
        'missing_pchembl', 'problematic_data_validity',
    }
    assert len(retained) == 0


def test_the_audit_keeps_the_excluded_rows():
    """What was dropped and why has to remain inspectable."""
    frame = pd.concat([_activities(activity_id=1),
                       _activities(activity_id=2, relationship_type='H')],
                      ignore_index=True)

    audit, retained = chembl.filter_activities(frame, direct_only=True)

    assert len(audit) == 2, 'the audit is the full input'
    assert len(retained) == 1
    assert sorted(audit.filter_status) == ['excluded', 'retained']


def test_filtering_does_not_modify_the_input_frame():
    frame = _activities(relationship_type='H')
    before = frame.copy()

    chembl.filter_activities(frame, direct_only=True)

    pd.testing.assert_frame_equal(frame, before)


# --------------------------------------------------------- evidence labelling

@pytest.mark.parametrize('origin,basis,exact,expected', [
    ('direct', 'canonical_sequence', True, 'direct_identical_sequence'),
    ('direct', 'canonical_sequence', False, 'direct_similar_sequence'),
    ('homology', 'variant_sequence', True, 'homology_identical_sequence'),
    ('homology', 'variant_sequence', False, 'homology_similar_sequence'),
])
def test_the_evidence_label_states_origin_and_sequence_agreement(origin, basis, exact, expected):
    row = {'assignment_origin': origin, 'sequence_basis': basis, 'exact_sequence_match': exact}

    assert chembl.evidence_level(row) == expected


def test_an_unresolved_variant_sequence_is_said_so_rather_than_called_similar():
    """
    The measurement was made on a mutant whose sequence ChEMBL does not provide,
    so neither 'identical' nor 'similar' would be true.
    """
    row = {'assignment_origin': 'direct',
           'sequence_basis': 'canonical_reference_variant_sequence_unavailable',
           'exact_sequence_match': True}

    assert chembl.evidence_level(row) == 'direct_variant_sequence_unresolved'


# ------------------------------------------------------------ donor mapping

def _targets_frame(rows):
    return pd.DataFrame(rows, columns=['tid', 'chembl_id', 'tax_id', 'target_type', 'pref_name'])


def _components_frame(rows):
    return pd.DataFrame(rows, columns=['tid', 'component_id', 'accession', 'sequence',
                                       'component_type'])


def _empty_activities():
    return pd.DataFrame(columns=['tid', 'variant_id', 'variant_mutation', 'variant_sequence'])


def test_a_single_protein_component_becomes_one_donor():
    targets = _targets_frame([(1, 'CHEMBL1', 2104, 'SINGLE PROTEIN', 'DHFR')])
    components = _components_frame([(1, 10, 'P12345', 'mka ilv*', 'PROTEIN')])

    donors = chembl.mapping_donors(targets, components, _empty_activities())

    assert len(donors) == 1
    assert donors.iloc[0]['sequence'] == 'MKAILV'
    assert donors.iloc[0]['sequence_basis'] == 'canonical_sequence'
    assert donors.iloc[0]['accession'] == 'P12345'


def test_a_target_with_several_protein_components_is_skipped(caplog):
    """
    A complex has no single sequence to assign, so guessing one would attach a
    measurement to the wrong protein.
    """
    targets = _targets_frame([(1, 'CHEMBL1', 2104, 'PROTEIN COMPLEX', 'complex')])
    components = _components_frame([
        (1, 10, 'P1', 'MKAILV', 'PROTEIN'),
        (1, 11, 'P2', 'MSNQVL', 'PROTEIN'),
    ])

    with caplog.at_level('WARNING'):
        donors = chembl.mapping_donors(targets, components, _empty_activities())

    assert donors.empty
    assert any('protein components' in record.message for record in caplog.records)


def test_a_target_with_no_protein_component_is_skipped():
    targets = _targets_frame([(1, 'CHEMBL1', 2104, 'SINGLE PROTEIN', 'x')])
    components = _components_frame([(1, 10, 'N1', 'ACGT', 'NUCLEIC ACID')])

    assert chembl.mapping_donors(targets, components, _empty_activities()).empty


def test_a_variant_adds_a_donor_with_its_own_sequence():
    targets = _targets_frame([(1, 'CHEMBL1', 2104, 'SINGLE PROTEIN', 'DHFR')])
    components = _components_frame([(1, 10, 'P12345', 'MKAILV', 'PROTEIN')])
    activities = pd.DataFrame([
        {'tid': 1, 'variant_id': 7, 'variant_mutation': 'A10V', 'variant_sequence': 'MKVILV'},
    ])

    donors = chembl.mapping_donors(targets, components, activities)

    assert len(donors) == 2, 'the canonical donor and the variant'
    variant = donors[donors.variant_id.notna()].iloc[0]
    assert variant['sequence'] == 'MKVILV'
    assert variant['sequence_basis'] == 'variant_sequence'
    assert variant['mutation'] == 'A10V'


def test_a_variant_without_a_sequence_falls_back_and_says_so():
    targets = _targets_frame([(1, 'CHEMBL1', 2104, 'SINGLE PROTEIN', 'DHFR')])
    components = _components_frame([(1, 10, 'P12345', 'MKAILV', 'PROTEIN')])
    activities = pd.DataFrame([
        {'tid': 1, 'variant_id': 7, 'variant_mutation': 'A10V', 'variant_sequence': None},
    ])

    donors = chembl.mapping_donors(targets, components, activities)

    variant = donors[donors.variant_id.notna()].iloc[0]
    assert variant['sequence'] == 'MKAILV', 'the canonical sequence is used'
    assert variant['sequence_basis'] == 'canonical_reference_variant_sequence_unavailable'


# --------------------------------------------------------- mapping guards

def _fasta(tmp_path, records):
    path = tmp_path / 'proteins.faa'
    path.write_text(''.join(f'>{name}\n{sequence}\n' for name, sequence in records))
    return path


@pytest.mark.parametrize('identity,coverage,threads', [
    (0, 70, 4), (101, 70, 4), (90, 0, 4), (90, 101, 4), (90, 70, 0), (90, 70, -1),
])
def test_out_of_range_mapping_parameters_are_refused(tmp_path, identity, coverage, threads):
    fasta = _fasta(tmp_path, [('GENE1', 'MKAILV')])
    donors = pd.DataFrame(columns=chembl.DONOR_COLUMNS)

    with pytest.raises(ValueError, match='Identity/coverage'):
        chembl.map_proteins(fasta, donors, tmp_path / 'work',
                            min_identity=identity, min_coverage=coverage, threads=threads)


def test_duplicate_locus_tags_are_refused(tmp_path):
    """A duplicated identifier would silently merge two genes' evidence."""
    fasta = _fasta(tmp_path, [('GENE1', 'MKAILV'), ('GENE1', 'MSNQVL')])
    donors = pd.DataFrame(columns=chembl.DONOR_COLUMNS)

    with pytest.raises(ValueError, match='unique locus_tag'):
        chembl.map_proteins(fasta, donors, tmp_path / 'work')


def test_an_empty_fasta_is_refused(tmp_path):
    fasta = tmp_path / 'proteins.faa'
    fasta.write_text('')
    donors = pd.DataFrame(columns=chembl.DONOR_COLUMNS)

    with pytest.raises(ValueError, match='unique locus_tag'):
        chembl.map_proteins(fasta, donors, tmp_path / 'work')


# ------------------------------------------------------------------ ligands

def test_compounds_are_aggregated_per_parent_molecule():
    activities = pd.DataFrame([
        {'parent_molregno': 1, 'compound_id': 'CHEMBL25', 'compound_name': 'aspirin',
         'canonical_smiles': 'CC(=O)O', 'standard_inchi_key': 'KEY1',
         'original_compound_id': 'CHEMBL25', 'tid': 1, 'activity_id': 100,
         'assignment_origin': 'direct'},
        {'parent_molregno': 1, 'compound_id': 'CHEMBL25', 'compound_name': 'aspirin',
         'canonical_smiles': 'CC(=O)O', 'standard_inchi_key': 'KEY1',
         'original_compound_id': 'CHEMBL26', 'tid': 2, 'activity_id': 101,
         'assignment_origin': 'homology'},
    ])

    ligands = chembl._ligands(activities)

    assert len(ligands) == 1
    row = ligands.iloc[0]
    assert row['bacterial_target_count'] == 2
    assert row['bacterial_measurement_count'] == 2
    assert row['bacterial_direct_measurement_count'] == 1
    assert row['bacterial_homology_measurement_count'] == 1
    assert row['original_compound_ids'] == 'CHEMBL25;CHEMBL26'


def test_rows_without_a_compound_are_not_ligands():
    activities = pd.DataFrame([
        {'parent_molregno': 1, 'compound_id': None, 'compound_name': None,
         'canonical_smiles': None, 'standard_inchi_key': None,
         'original_compound_id': None, 'tid': 1, 'activity_id': 100,
         'assignment_origin': 'direct'},
    ])

    assert chembl._ligands(activities).empty
