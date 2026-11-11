"""The UniParc exact-identity tier, with the HTTP responses supplied locally.

This tier exists to recover genes whose UniProt entry was pruned while AlphaFold
DB still serves a model of the identical sequence. Its whole justification is
that it cannot return a homolog: UniParc is keyed on the sequence, and the model
is accepted only when its sequence is byte-identical to the query. These tests
are mostly about the refusals, because an accepted wrong mapping would silently
attach another protein's structure to a gene.

No network: every request is answered from a table keyed by URL.
"""

import json
import os

import pytest

from ftscripts import structures


SPECIES, STRAIN = 2104, 272634
SEQUENCE = 'MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCK'
CHECKSUM = structures.protein_sequence_checksum(SEQUENCE)
UPI = 'UPI0000000001'


class FakeResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def json(self):
        if isinstance(self._payload, Exception):
            raise self._payload
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            import requests
            raise requests.HTTPError(f'status {self.status_code}')


@pytest.fixture
def http(monkeypatch):
    """
    Answer requests.get from a routing table, and record what was asked for.

    Routes are matched on a substring of the URL. A route may be a payload, an
    (payload, status) pair, or an exception to raise.
    """
    state = {'routes': {}, 'calls': []}

    def fake_get(url, params=None, timeout=None, headers=None):
        state['calls'].append((url, params))
        for fragment, route in state['routes'].items():
            if fragment in url:
                if isinstance(route, BaseException):
                    raise route
                if isinstance(route, tuple):
                    payload, status = route
                    return FakeResponse(payload, status)
                return FakeResponse(route)
        return FakeResponse(None, 404)

    monkeypatch.setattr(structures.requests, 'get', fake_get)
    # No real waiting between retries.
    monkeypatch.setattr(structures.time, 'sleep', lambda seconds: None)
    return state


def _search_hit(upi=UPI):
    return {'results': [{'uniParcId': upi}]}


def _entry(xrefs):
    return {'uniParcId': UPI, 'uniParcCrossReferences': xrefs}


def _xref(accession, taxid, active=True):
    return {'database': 'UniProtKB/TrEMBL', 'id': accession, 'active': active,
            'organism': {'taxonId': taxid}}


def _model(sequence, accession='A00001', plddt=78.4, version=6):
    return [{
        'uniprotSequence': sequence,
        'globalMetricValue': plddt,
        'latestVersion': version,
        'cifUrl': f'https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.cif',
        'pdbUrl': f'https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v{version}.pdb',
    }]


# ----------------------------------------------------------- the HTTP helper

def test_absent_resources_are_not_retried(http):
    """404, 400 and 410 are answers, not failures: retrying them wastes time."""
    for status in (400, 404, 410):
        http['routes'] = {'example': (None, status)}
        http['calls'].clear()

        assert structures._rest_get_json('https://example.org/thing', retries=3) is None
        assert len(http['calls']) == 1, f'status {status} was retried'


def test_a_transient_failure_is_retried_and_can_succeed(http, monkeypatch):
    import requests

    attempts = {'n': 0}

    def flaky(url, params=None, timeout=None, headers=None):
        attempts['n'] += 1
        if attempts['n'] < 3:
            raise requests.ConnectionError('reset')
        return FakeResponse({'ok': True})

    monkeypatch.setattr(structures.requests, 'get', flaky)

    assert structures._rest_get_json('https://example.org/x', retries=3) == {'ok': True}
    assert attempts['n'] == 3


def test_giving_up_returns_none_rather_than_raising(http):
    import requests

    http['routes'] = {'example': requests.ConnectionError('down')}

    assert structures._rest_get_json('https://example.org/x', retries=2) is None


def test_a_server_error_is_retried_then_gives_up(http):
    http['routes'] = {'example': (None, 500)}

    assert structures._rest_get_json('https://example.org/x', retries=2) is None
    assert len(http['calls']) == 2


def test_an_unparseable_body_is_not_a_crash(http):
    http['routes'] = {'example': ValueError('not json')}

    assert structures._rest_get_json('https://example.org/x', retries=1) is None


# --------------------------------------------------- checksum to UniParc entry

def test_a_checksum_resolves_through_the_entry_endpoint(http):
    """
    The search endpoint omits cross-references, so the entry endpoint has to be
    fetched: ranking by strain needs the organism and the active flag per xref.
    """
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([_xref('A00001', STRAIN)]),
    }

    entry = structures.uniparc_entry_by_checksum(CHECKSUM)

    assert entry['uniParcId'] == UPI
    assert entry['uniParcCrossReferences'][0]['id'] == 'A00001'
    assert any('checksum:' in str(params) for _, params in http['calls'] if params)


def test_a_sequence_absent_from_uniparc_maps_to_nothing(http):
    http['routes'] = {'/uniparc/search': {'results': []}}

    assert structures.uniparc_entry_by_checksum(CHECKSUM) is None


def test_the_search_hit_is_used_when_the_entry_has_no_cross_references(http):
    http['routes'] = {
        '/uniparc/search': {'results': [{'uniParcId': UPI, 'sequence': {'value': SEQUENCE}}]},
        f'/uniparc/{UPI}.json': None,
    }

    entry = structures.uniparc_entry_by_checksum(CHECKSUM)

    assert entry['uniParcId'] == UPI


# ---------------------------------------------------------------- AFDB models

def test_a_model_is_reported_with_its_sequence_and_confidence(http):
    http['routes'] = {'/api/prediction/A00001': _model(SEQUENCE, 'A00001', plddt=91.2)}

    model = structures.alphafold_model_for_accession('A00001')

    assert model['sequence'] == SEQUENCE
    assert model['plddt'] == 91.2
    assert model['version'] == 6
    assert model['pdb_url'].endswith('.pdb')


def test_an_accession_without_a_model_is_none(http):
    http['routes'] = {'/api/prediction/A00001': (None, 404)}

    assert structures.alphafold_model_for_accession('A00001') is None


def test_an_unexpected_payload_shape_is_none(http):
    """The API returns a list; anything else is not a model."""
    http['routes'] = {'/api/prediction/A00001': {'error': 'nope'}}

    assert structures.alphafold_model_for_accession('A00001') is None


# -------------------------------------------------- the mapping decision

def test_an_identical_model_from_the_strain_is_accepted(http):
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([_xref('STRAINACC', STRAIN)]),
        '/api/prediction/STRAINACC': _model(SEQUENCE, 'STRAINACC'),
    }

    mapping = structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN)

    assert mapping['uniprot_id'] == 'STRAINACC'
    assert mapping['uniparc_id'] == UPI
    assert mapping['checksum'] == CHECKSUM
    assert mapping['tier'] == 4


def test_a_model_whose_sequence_differs_is_refused(http):
    """
    The point of the tier. A model is built for one sequence version, so if the
    entry changed after the model was made, the model is of another protein.
    """
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([_xref('STRAINACC', STRAIN)]),
        '/api/prediction/STRAINACC': _model(SEQUENCE[:-1] + 'G', 'STRAINACC'),
    }

    assert structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN) is None


def test_a_model_one_residue_shorter_is_refused(http):
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([_xref('STRAINACC', STRAIN)]),
        '/api/prediction/STRAINACC': _model(SEQUENCE[:-1], 'STRAINACC'),
    }

    assert structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN) is None


def test_formatting_differences_do_not_cause_a_refusal(http):
    """Identity is of residues, not of whitespace or a trailing stop."""
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([_xref('STRAINACC', STRAIN)]),
        '/api/prediction/STRAINACC': _model(f'{SEQUENCE.lower()}*', 'STRAINACC'),
    }

    mapping = structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN)

    assert mapping['uniprot_id'] == 'STRAINACC'


def test_candidates_are_tried_in_tier_order_and_the_first_match_wins(http):
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([
            _xref('OTHERACC', 9606),
            _xref('SPECIESACC', SPECIES),
            _xref('STRAINACC', STRAIN),
        ]),
        '/api/prediction/STRAINACC': _model(SEQUENCE, 'STRAINACC'),
        '/api/prediction/SPECIESACC': _model(SEQUENCE, 'SPECIESACC'),
        '/api/prediction/OTHERACC': _model(SEQUENCE, 'OTHERACC'),
    }

    mapping = structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN)

    assert mapping['uniprot_id'] == 'STRAINACC', 'the strain must be preferred'
    asked = [url for url, _ in http['calls'] if 'prediction' in url]
    assert asked[0].endswith('STRAINACC'), 'the strain was not tried first'


def test_a_lower_tier_is_used_when_the_closer_one_has_no_model(http):
    """
    This is the realistic case: the strain accession was pruned and dropped from
    AFDB, while a species-level accession for the identical sequence still has a
    model.
    """
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([
            _xref('STRAINACC', STRAIN, active=False),
            _xref('SPECIESACC', SPECIES),
        ]),
        '/api/prediction/STRAINACC': (None, 404),
        '/api/prediction/SPECIESACC': _model(SEQUENCE, 'SPECIESACC'),
    }

    mapping = structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN)

    assert mapping['uniprot_id'] == 'SPECIESACC'
    assert mapping['tier'] == 5


def test_no_candidate_with_a_model_maps_to_nothing(http):
    """
    An old deletion loses the models too, which is why ColabFold is the only
    tier with a floor.
    """
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([
            _xref('GONE1', STRAIN, active=False),
            _xref('GONE2', SPECIES, active=False),
        ]),
        '/api/prediction/GONE1': (None, 404),
        '/api/prediction/GONE2': (None, 404),
    }

    assert structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN) is None


def test_an_empty_sequence_is_not_looked_up(http):
    http['routes'] = {}

    assert structures.uniparc_mapping_for_sequence('   ', SPECIES, STRAIN) is None
    assert http['calls'] == [], 'nothing should be requested for an empty sequence'


def test_the_accepted_mapping_carries_the_model_metadata(http):
    http['routes'] = {
        '/uniparc/search': _search_hit(),
        f'/uniparc/{UPI}.json': _entry([_xref('STRAINACC', STRAIN)]),
        '/api/prediction/STRAINACC': _model(SEQUENCE, 'STRAINACC', plddt=64.5, version=6),
    }

    mapping = structures.uniparc_mapping_for_sequence(SEQUENCE, SPECIES, STRAIN)

    assert mapping['plddt'] == 64.5
    assert mapping['alphafold_version'] == 6
    assert mapping['active'] is True
    assert mapping['taxid'] == str(STRAIN)
