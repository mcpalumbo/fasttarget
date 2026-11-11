"""The fallback as a pipeline stage: per-gene recovery and its annotations.

The tier runs after BLAST, so its contract is narrow: never displace a mapping
that already succeeded, never take down the run when a lookup fails, and produce
annotation records for accessions UniProtKB will no longer describe.
"""

import json
import os

import pytest

from conftest import GENES, SEQUENCES
from ftscripts import structures


SPECIES, STRAIN = 2104, 272634


@pytest.fixture
def mapping_calls(monkeypatch):
    """Replace the per-sequence lookup, recording what it was asked about."""
    state = {'asked': [], 'results': {}}

    def fake_mapping(sequence, specie_taxid, strain_taxid, timeout=30, retries=3):
        state['asked'].append(sequence)
        return state['results'].get(structures.normalize_protein_sequence(sequence))

    monkeypatch.setattr(structures, 'uniparc_mapping_for_sequence', fake_mapping)
    return state


def _result(accession, tier=4, plddt=80.0, taxid=STRAIN, sequence=None):
    return {
        'uniprot_id': accession,
        'uniparc_id': 'UPI0000000001',
        'checksum': 'AAAAAAAAAAAAAAAA',
        'tier': tier,
        'taxid': str(taxid),
        'active': False,
        'plddt': plddt,
        'alphafold_version': 6,
        'sequence': sequence,
    }


def _fallback_path(output_path, organism_name):
    return os.path.join(output_path, organism_name, 'structures', 'uniprot_files',
                        f'uniparc_{organism_name}_fallback.json')


# ------------------------------------------------------------------- mapping

def test_only_the_unmapped_genes_are_looked_up(organism, mapping_calls):
    """The tier runs after BLAST and must not touch genes already mapped."""
    output_path, organism_name = organism
    unmapped = [GENES[1], GENES[3]]
    mapping_calls['results'] = {
        structures.normalize_protein_sequence(SEQUENCES[gene]): _result(f'ACC_{gene}')
        for gene in unmapped
    }

    recovered = structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, unmapped
    )

    assert sorted(recovered) == sorted(unmapped)
    assert len(mapping_calls['asked']) == len(unmapped)
    for sequence in mapping_calls['asked']:
        assert sequence in SEQUENCES.values()


def test_no_unmapped_genes_means_no_work(organism, mapping_calls):
    output_path, organism_name = organism

    assert structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, []
    ) == {}
    assert mapping_calls['asked'] == []


def test_genes_that_resolve_to_nothing_are_simply_absent(organism, mapping_calls):
    output_path, organism_name = organism
    unmapped = [GENES[0], GENES[1]]
    mapping_calls['results'] = {
        structures.normalize_protein_sequence(SEQUENCES[GENES[0]]): _result('ACC_0')
    }

    recovered = structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, unmapped
    )

    assert list(recovered) == [GENES[0]]


def test_a_failing_lookup_does_not_abort_the_other_genes(organism, monkeypatch):
    """A fallback tier must never take down a mapping that already succeeded."""
    output_path, organism_name = organism

    def sometimes_raises(sequence, specie_taxid, strain_taxid, timeout=30, retries=3):
        if sequence == SEQUENCES[GENES[0]]:
            raise RuntimeError('the API misbehaved')
        return _result('ACC_SURVIVOR')

    monkeypatch.setattr(structures, 'uniparc_mapping_for_sequence', sometimes_raises)

    recovered = structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, [GENES[0], GENES[1]]
    )

    assert list(recovered) == [GENES[1]]


def test_a_gene_absent_from_the_proteome_is_skipped(organism, mapping_calls):
    output_path, organism_name = organism

    recovered = structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, ['NOT_A_GENE']
    )

    assert recovered == {}
    assert mapping_calls['asked'] == []


def test_a_missing_proteome_skips_the_tier_without_failing(organism, mapping_calls):
    """
    The tier is optional, so a missing FASTA is a reason to skip it rather than
    to stop a run whose other tiers worked.
    """
    output_path, organism_name = organism
    os.remove(os.path.join(output_path, organism_name, 'genome', f'{organism_name}.faa'))

    recovered = structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, [GENES[0]]
    )

    assert recovered == {}
    assert mapping_calls['asked'] == []


def test_the_recovered_mapping_is_written_for_later_stages(organism, mapping_calls):
    output_path, organism_name = organism
    mapping_calls['results'] = {
        structures.normalize_protein_sequence(SEQUENCES[GENES[0]]): _result('ACC_0')
    }

    structures.uniparc_fallback_mapping(
        output_path, organism_name, SPECIES, STRAIN, [GENES[0]]
    )

    with open(_fallback_path(output_path, organism_name)) as handle:
        written = json.load(handle)
    assert written[GENES[0]]['uniprot_id'] == 'ACC_0'


# --------------------------------------------------------------- annotations

def _write_fallback(output_path, organism_name, payload):
    path = _fallback_path(output_path, organism_name)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        json.dump(payload, handle)
    return path


def test_annotations_are_synthesised_for_withdrawn_accessions(organism):
    """
    UniProtKB serves a withdrawn accession as an empty document, so without this
    the gene would be reported as having no structure at all.
    """
    output_path, organism_name = organism
    sequence = SEQUENCES[GENES[0]]
    _write_fallback(output_path, organism_name, {
        GENES[0]: _result('A0AB38W5R8', tier=5, plddt=78.38, taxid=SPECIES,
                          sequence=sequence),
    })

    annotations = structures.uniparc_fallback_annotations(output_path, organism_name)

    record = annotations['A0AB38W5R8']
    assert record['AlphaFoldDB'] == 'A0AB38W5R8'
    assert record['Locus_tag'] == GENES[0]
    assert record['Sequence'] == sequence
    assert record['PDB_id'] is None, 'a predicted model is not an experimental structure'
    assert record['Dataset'] == 'UniParc fallback'
    assert record['UniParc_tier'] == 5
    assert record['AlphaFold_plddt'] == 78.38


def test_the_provenance_of_each_recovered_gene_is_recorded(organism):
    output_path, organism_name = organism
    _write_fallback(output_path, organism_name, {
        GENES[0]: _result('ACC_A', tier=4, taxid=STRAIN, sequence=SEQUENCES[GENES[0]]),
        GENES[1]: _result('ACC_B', tier=5, taxid=SPECIES, sequence=SEQUENCES[GENES[1]]),
    })

    annotations = structures.uniparc_fallback_annotations(output_path, organism_name)

    assert annotations['ACC_A']['UniParc_taxid'] == str(STRAIN)
    assert annotations['ACC_B']['UniParc_taxid'] == str(SPECIES)
    assert annotations['ACC_A']['UniParc_id'] == 'UPI0000000001'
    assert annotations['ACC_A']['UniProt_active'] is False


def test_no_fallback_file_means_no_annotations(organism):
    output_path, organism_name = organism

    assert structures.uniparc_fallback_annotations(output_path, organism_name) == {}


def test_an_empty_fallback_file_means_no_annotations(organism):
    output_path, organism_name = organism
    _write_fallback(output_path, organism_name, {})

    assert structures.uniparc_fallback_annotations(output_path, organism_name) == {}
