"""AlphaFill: accession safety, confidence tiers and metadata parsing.

The confidence tiers decide which transplanted ligands are reported as usable,
so their boundaries are the thing to pin: a limit that is inclusive on the wrong
side moves ligands between tiers silently.
"""

import json
import math
import os
from pathlib import Path

import pandas as pd
import pytest

from ftscripts import alphafill


# --------------------------------------------------- identifiers and paths

@pytest.mark.parametrize('value,expected', [
    ('P12345', 'P12345'),
    ('p12345', 'P12345'),
    ('  P12345  ', 'P12345'),
    ('A0AB38W5R8', 'A0AB38W5R8'),
    ('A0A0K2_X-1', 'A0A0K2_X-1'),
])
def test_an_accession_is_normalised_to_upper_case(value, expected):
    assert alphafill.sanitize_uniprot_id(value) == expected


@pytest.mark.parametrize('value', [
    '', '   ', None, '..', '../../etc/passwd', 'P12345/../x', 'P/12345',
    'P12345\0', 'x' * 33,
])
def test_an_accession_that_could_escape_its_directory_is_refused(value):
    """The accession becomes both a URL path and a directory name."""
    with pytest.raises(ValueError):
        alphafill.sanitize_uniprot_id(value)


@pytest.mark.parametrize('value', ['MPN_RS00140', 'SYN_0001', 'gene.1'])
def test_a_locus_tag_keeps_its_own_characters(value):
    assert alphafill._validate_locus_tag(value) == value


@pytest.mark.parametrize('value', ['', '  ', '.', '..', 'a/b', 'a\\b', 'a\0b', None, 42])
def test_a_locus_tag_that_is_not_a_single_directory_name_is_refused(value):
    with pytest.raises(ValueError):
        alphafill._validate_locus_tag(value)


def test_artifact_paths_are_stored_relative_to_the_alphafill_root(tmp_path):
    """
    Artifact paths are recorded relative to ligands/alphafill/, so an output
    directory that moves does not invalidate them.
    """
    root = Path(str(tmp_path)) / 'ligands' / 'alphafill'
    root.mkdir(parents=True)
    inside = root / 'models' / 'GENE' / 'P1.cif'

    converted = alphafill._portable_paths({'cif_path': str(inside)}, root)

    assert converted['cif_path'] == 'models/GENE/P1.cif'


def test_only_the_named_artifact_fields_are_converted(tmp_path):
    """URLs and diagnostics are not local paths and must survive untouched."""
    root = Path(str(tmp_path)) / 'ligands' / 'alphafill'
    root.mkdir(parents=True)
    url = 'https://alphafill.eu/v1/aff/P12345'

    converted = alphafill._portable_paths(
        {'cif_path': str(root / 'a.cif'), 'source_url': url, 'error_message': '/tmp/x'}, root
    )

    assert converted['cif_path'] == 'a.cif'
    assert converted['source_url'] == url
    assert converted['error_message'] == '/tmp/x'


def test_nested_structures_are_converted_throughout(tmp_path):
    root = Path(str(tmp_path)) / 'ligands' / 'alphafill'
    root.mkdir(parents=True)

    converted = alphafill._portable_paths(
        {'results': [{'metadata_path': str(root / 'GENE' / 'm.json')}]}, root
    )

    assert converted['results'][0]['metadata_path'] == 'GENE/m.json'


# ------------------------------------------------------------------- numbers

@pytest.mark.parametrize('value,expected', [
    (1, 1.0), ('2.5', 2.5), (0, 0.0), (-3, -3.0),
    (None, None), (True, None), (False, None),
    ('abc', None), ('', None), (float('nan'), None),
    (float('inf'), None), (float('-inf'), None),
])
def test_only_finite_numbers_are_accepted(value, expected):
    """A boolean is not a measurement, and neither is NaN or infinity."""
    assert alphafill._number(value) == expected


@pytest.mark.parametrize('value,expected', [
    (0.5, 50.0),      # a fraction is scaled to percent
    (0.0, 0.0),
    (1.0, 100.0),     # the ambiguous point: treated as a fraction
    (1.5, 1.5),       # above one, already a percentage
    (45.0, 45.0),
    (100.0, 100.0),
    (-1, None),       # outside 0-100
    (101, None),
    (None, None),
    ('abc', None),
])
def test_identity_is_normalised_to_percent(value, expected):
    assert alphafill.normalize_identity(value) == expected


def test_an_identity_of_exactly_one_is_read_as_a_fraction():
    """
    Documenting the ambiguity rather than hiding it: AlphaFill reports identity
    either as a fraction or a percentage, so 1 means 100% and a genuine 1%
    identity cannot be expressed. Such a donor would be rejected on structural
    metrics anyway.
    """
    assert alphafill.normalize_identity(1) == 100.0


# -------------------------------------------------------- identity binning

@pytest.mark.parametrize('identity,expected', [
    (0, 'below_25'), (24.9, 'below_25'),
    (25, '25-35'), (34.9, '25-35'),
    (35, '35-45'), (44.9, '35-45'),
    (45, '45-55'), (55, '55-65'), (65, '65-75'), (74.9, '65-75'),
    (75, '>=75'), (100, '>=75'),
    (None, 'unknown'), (-1, 'unknown'), (101, 'unknown'), ('abc', 'unknown'),
])
def test_identity_bins_are_closed_on_the_lower_bound(identity, expected):
    assert alphafill.identity_category(identity) == expected


def test_every_identity_bin_is_reachable():
    """A bin no value can fall into would be a silent dead category."""
    produced = {alphafill.identity_category(value)
                for value in [10, 30, 40, 50, 60, 70, 80, None]}
    assert produced == set(alphafill.IDENTITY_CATEGORIES)


# ------------------------------------------------------ confidence tiers

def _ligand(local_rmsd, tcs, identity=50.0):
    return {'compound_code': 'ATP', 'local_rmsd': local_rmsd, 'tcs': tcs,
            'donor_identity': identity}


@pytest.mark.parametrize('local_rmsd,tcs,expected', [
    (0.5, 0.5, 'high'),
    (0.91, 0.63, 'high'),
    (0.92, 0.5, 'intermediate'),   # the limit is strict: at 0.92 it is not high
    (0.5, 0.64, 'intermediate'),
    (3.09, 1.26, 'intermediate'),
    (3.10, 0.5, 'low'),            # at the intermediate limit it is low
    (0.5, 1.27, 'low'),
    (10.0, 5.0, 'low'),
])
def test_both_metrics_must_satisfy_a_tier(local_rmsd, tcs, expected):
    assert alphafill.filter_ligand(_ligand(local_rmsd, tcs))['confidence_category'] == expected


@pytest.mark.parametrize('local_rmsd,tcs', [
    (None, 0.5), (0.5, None), (None, None), ('abc', 0.5), (float('nan'), 0.5),
])
def test_a_missing_metric_means_unknown_not_low(local_rmsd, tcs):
    """Unknown and low are different: one was not measured, the other failed."""
    result = alphafill.filter_ligand(_ligand(local_rmsd, tcs))

    assert result['confidence_category'] == 'unknown'
    assert any(reason.startswith('missing_') for reason in result['rejection_reasons'])


def test_a_negative_metric_is_invalid_rather_than_excellent():
    """Without this check a negative RMSD would pass every threshold."""
    result = alphafill.filter_ligand(_ligand(-1.0, -1.0))

    assert result['confidence_category'] == 'unknown'
    assert 'invalid_local_rmsd' in result['rejection_reasons']
    assert 'invalid_tcs' in result['rejection_reasons']


@pytest.mark.parametrize('confidence_input,accepted', [
    ((0.5, 0.5), True), ((3.0, 1.0), True), ((5.0, 5.0), False), ((None, None), False),
])
def test_only_high_and_intermediate_are_accepted(confidence_input, accepted):
    local_rmsd, tcs = confidence_input
    result = alphafill.filter_ligand(_ligand(local_rmsd, tcs))

    assert (result['filter_status'] == 'accepted') is accepted


def test_identity_and_global_rmsd_never_reject_a_transplant():
    """They are reported for judgement, not used as filters."""
    ligand = dict(_ligand(0.5, 0.5, identity=5.0), global_rmsd=99.0)

    result = alphafill.filter_ligand(ligand)

    assert result['confidence_category'] == 'high'
    assert result['identity_category'] == 'below_25'


def test_the_thresholds_can_be_tightened(tmp_path):
    ligand = _ligand(0.8, 0.5)

    assert alphafill.filter_ligand(ligand)['confidence_category'] == 'high'
    assert alphafill.filter_ligand(
        ligand, high_local_rmsd=0.5, intermediate_local_rmsd=1.0
    )['confidence_category'] == 'intermediate'


# ---------------------------------------------------- compound summaries

def test_a_compound_is_counted_once_at_its_best_confidence():
    """
    Donor structure, chain and instance do not define a new compound: the same
    ligand placed from three donors is one compound.
    """
    ligands = [
        {'compound_code': 'ATP', 'confidence_category': 'low'},
        {'compound_code': 'ATP', 'confidence_category': 'high'},
        {'compound_code': 'ATP', 'confidence_category': 'intermediate'},
        {'compound_code': 'GTP', 'confidence_category': 'intermediate'},
    ]

    summary = alphafill.summarize_compounds(ligands)

    assert summary['alphafill_high_confidence_compounds'] == 'ATP'
    assert summary['alphafill_intermediate_confidence_compounds'] == 'GTP'
    assert summary['alphafill_low_confidence_compounds'] == ''
    assert summary['alphafill_unique_compounds'] == 2
    assert summary['alphafill_unique_accepted_compounds'] == 2


def test_compound_codes_are_normalised_and_blanks_ignored():
    ligands = [
        {'compound_code': ' atp '}, {'compound_code': 'ATP'},
        {'compound_code': ''}, {'compound_code': None}, {},
    ]

    summary = alphafill.summarize_compounds(ligands)

    assert summary['alphafill_unique_compounds'] == 1
    assert summary['alphafill_unknown_confidence_compounds'] == 'ATP'


def test_an_unusable_result_has_null_counts_not_zero():
    """
    A failed search and a successful search with no hits must not look alike:
    zero is a finding, null is the absence of one.
    """
    empty_but_usable = alphafill.summarize_compounds([], usable=True)
    unusable = alphafill.summarize_compounds([], usable=False)

    assert empty_but_usable['alphafill_unique_compounds'] == 0
    assert unusable['alphafill_unique_compounds'] is None
    assert unusable['alphafill_high_confidence_compound_count'] is None


def test_an_unrecognised_confidence_is_treated_as_unknown():
    summary = alphafill.summarize_compounds([
        {'compound_code': 'ATP', 'confidence_category': 'excellent'},
    ])

    assert summary['alphafill_unknown_confidence_compounds'] == 'ATP'


# --------------------------------------------------------- metadata parsing

def _payload(transplants, pdb_id='1ABC', identity=0.45):
    return {'hits': [{
        'pdb_id': pdb_id,
        'pdb_asym_id': 'A',
        'alignment': {'identity': identity},
        'global_rmsd': 1.2,
        'transplants': transplants,
    }]}


def _transplant(analogue='ATP', asym='B', local_rmsd=0.5, tcs=0.3, **extra):
    return {'analogue_id': analogue, 'compound_id': analogue, 'asym_id': asym,
            'local_rmsd': local_rmsd, 'clash': {'score': tcs}, **extra}


def test_a_transplant_becomes_one_candidate_with_its_metrics():
    candidates = alphafill.parse_alphafill_metadata('p12345', _payload([_transplant()]))

    assert len(candidates) == 1
    candidate = candidates[0]
    assert candidate['uniprot_id'] == 'P12345'
    assert candidate['compound_code'] == 'ATP'
    assert candidate['donor_pdb'] == '1ABC'
    assert candidate['local_rmsd'] == 0.5
    assert candidate['tcs'] == 0.3
    assert candidate['donor_identity'] == 45.0


def test_no_hits_is_an_empty_result_not_an_error():
    assert alphafill.parse_alphafill_metadata('P12345', {'hits': None}) == []
    assert alphafill.parse_alphafill_metadata('P12345', {'hits': []}) == []


def test_a_hit_without_transplants_contributes_nothing():
    payload = {'hits': [{'pdb_id': '1ABC', 'transplants': None}]}

    assert alphafill.parse_alphafill_metadata('P12345', payload) == []


@pytest.mark.parametrize('payload,fragment', [
    ('not a dict', 'JSON object'),
    ({}, 'hits'),
    ({'hits': 'nope'}, 'list or null'),
    ({'hits': ['nope']}, 'object'),
    ({'hits': [{'pdb_id': '1ABC'}]}, 'transplants'),
    ({'hits': [{'transplants': 'nope'}]}, 'list or null'),
    ({'hits': [{'transplants': ['nope']}]}, 'object'),
])
def test_a_malformed_payload_is_refused_with_the_reason(payload, fragment):
    """
    Parsing must not guess: a changed API schema has to stop the stage rather
    than quietly produce fewer ligands.
    """
    with pytest.raises(ValueError, match=fragment):
        alphafill.parse_alphafill_metadata('P12345', payload)


@pytest.mark.parametrize('missing', ['analogue_id', 'compound_id', 'asym_id'])
def test_a_transplant_without_its_identifying_fields_is_refused(missing):
    transplant = _transplant()
    del transplant[missing]

    with pytest.raises(ValueError, match=missing):
        alphafill.parse_alphafill_metadata('P12345', _payload([transplant]))


def test_ligand_instance_ids_do_not_depend_on_array_order():
    """
    The IDs identify a placed ligand across runs, so they are derived from the
    instance and not from its position in the response.
    """
    first = _transplant(analogue='ATP', asym='B')
    second = _transplant(analogue='GTP', asym='C')

    one_way = alphafill.parse_alphafill_metadata('P12345', _payload([first, second]))
    other_way = alphafill.parse_alphafill_metadata('P12345', _payload([second, first]))

    assert [c['ligand_instance_id'] for c in one_way] == [c['ligand_instance_id'] for c in other_way]


def test_repeated_instances_get_distinct_ids():
    """Two copies of the same ligand in the same chain are two placements."""
    identical = _transplant(analogue='ATP', asym='B')

    candidates = alphafill.parse_alphafill_metadata(
        'P12345', _payload([dict(identical), dict(identical)])
    )

    ids = [candidate['ligand_instance_id'] for candidate in candidates]
    assert len(set(ids)) == 2


def test_the_placed_compound_is_preferred_over_the_donor_compound():
    """analogue_id is what was placed in the model; compound_id is the donor's."""
    transplant = _transplant()
    transplant['analogue_id'] = 'ADP'
    transplant['compound_id'] = 'ATP'

    candidate = alphafill.parse_alphafill_metadata('P12345', _payload([transplant]))[0]

    assert candidate['compound_code'] == 'ADP'
    assert candidate['donor_compound_id'] == 'ATP'


# --------------------------------------------------------------- the options

def test_the_documented_defaults_are_what_the_code_uses():
    options = alphafill.validate_options({})

    assert options['high_local_rmsd'] == 0.92
    assert options['high_tcs'] == 0.64
    assert options['intermediate_local_rmsd'] == 3.10
    assert options['intermediate_tcs'] == 1.27
    assert options['download_cif'] == 'accepted_only'
    assert options['retries'] == 3
    assert options['max_workers'] == 4


@pytest.mark.parametrize('key', ['min_identity', 'max_global_rmsd', 'max_local_rmsd'])
def test_a_retired_filter_key_is_refused_with_an_explanation(key):
    """Silently ignoring it would leave the user believing a filter applies."""
    with pytest.raises(ValueError, match='confidence categories'):
        alphafill.validate_options({key: 50})


def test_the_renamed_backoff_key_is_refused():
    with pytest.raises(ValueError, match='backoff'):
        alphafill.validate_options({'retry_backoff': 2})


def test_a_tier_limit_must_be_below_the_next_one():
    with pytest.raises(ValueError, match='must be less than'):
        alphafill.validate_options({'high_local_rmsd': 5.0, 'intermediate_local_rmsd': 3.0})


@pytest.mark.parametrize('config,fragment', [
    ({'high_tcs': 0}, 'positive finite'),
    ({'high_tcs': -1}, 'positive finite'),
    ({'timeout': 0}, 'timeout'),
    ({'retries': -1}, 'retries'),
    ({'backoff': -1}, 'backoff'),
    ({'max_workers': 0}, 'max_workers'),
    ({'max_workers': 65}, 'max_workers'),
    ({'download_cif': 'sometimes'}, 'download_cif'),
    ({'refresh_cache': 'yes'}, 'refresh_cache'),
    ({'strict_errors': 1}, 'strict_errors'),
])
def test_an_unusable_option_is_reported_by_name(config, fragment):
    with pytest.raises(ValueError, match=fragment):
        alphafill.validate_options(config)


def test_every_error_is_reported_not_only_the_first():
    with pytest.raises(ValueError) as error:
        alphafill.validate_options({'retries': -1, 'max_workers': 99, 'download_cif': 'maybe'})

    message = str(error.value)
    assert 'retries' in message and 'max_workers' in message and 'download_cif' in message
