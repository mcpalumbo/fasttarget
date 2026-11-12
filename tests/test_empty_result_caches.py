"""Empty results that are valid, and must not be cached forever.

A step with nothing to do writes '{}' or a zero-byte file. file_check reads two
bytes as content, so a run made before the inputs existed poisoned the cache: in
the twenty-protein case three Foldseek files written by a run without structures
kept every later run reporting no human structural hits, while 19 result files
with thousands of hit lines sat on disk unread.
"""

import json
import os

import pytest

from conftest import GENES
from ftscripts import files, offtargets


@pytest.fixture
def foldseek_results(organism):
    """A Foldseek run that produced real results for four of five genes."""
    output_path, organism_name = organism
    results_dir = os.path.join(output_path, organism_name, 'offtarget', 'foldseek_results')

    mapping = {}
    for gene in GENES[:4]:
        gene_dir = os.path.join(results_dir, f'AF_{gene}_output_foldseek')
        os.makedirs(gene_dir, exist_ok=True)
        path = os.path.join(gene_dir, f'AF_{gene}_vs_DB_human_reference_foldseek_results.tsv')
        with open(path, 'w') as handle:
            # The parser selects columns by name, so the header is part of the format.
            handle.write(
                'query\ttarget\talnlen\tqcov\ttcov\tlddt\tqtmscore\tttmscore\t'
                'alntmscore\trmsd\tprob\tpident\tevalue\n'
            )
            handle.write(
                f'AF_{gene}\tHUMAN_{gene}\t250\t0.95\t0.90\t0.70\t0.85\t0.80\t0.88\t2.1\t0.99\t34.5\t1e-40\n'
            )
        mapping[gene] = [path]
    mapping[GENES[4]] = []

    return output_path, organism_name, results_dir, mapping


def test_an_empty_parse_is_not_reused_once_results_appear(foldseek_results):
    output_path, organism_name, results_dir, mapping = foldseek_results

    # What a run with no structures leaves behind.
    dict_path = os.path.join(results_dir, 'human_foldseek_dict.json')
    with open(dict_path, 'w') as handle:
        json.dump({}, handle)
    assert files.file_check(dict_path) is True, 'two bytes look like content'

    parsed = offtargets.foldseek_human_parser(output_path, organism_name, mapping)

    assert len(parsed) == len(mapping), 'the stale empty dict was reused'
    with_hit = [gene for gene, value in parsed.items() if value.get('target_foldseek')]
    assert sorted(with_hit) == sorted(GENES[:4])


def test_an_empty_parse_is_reused_when_there_is_still_nothing_to_parse(organism):
    output_path, organism_name = organism
    results_dir = os.path.join(output_path, organism_name, 'offtarget', 'foldseek_results')
    os.makedirs(results_dir, exist_ok=True)

    first = offtargets.foldseek_human_parser(output_path, organism_name, {})
    assert first == {}

    # Second call with the same (empty) input must be a cache hit, not a rebuild.
    dict_path = os.path.join(results_dir, 'human_foldseek_dict.json')
    before = os.path.getmtime(dict_path)
    second = offtargets.foldseek_human_parser(output_path, organism_name, {})
    assert second == {}
    assert os.path.getmtime(dict_path) == before


def test_the_merge_step_is_not_reused_when_the_parse_changes(foldseek_results):
    output_path, organism_name, results_dir, mapping = foldseek_results
    offtarget_dir = os.path.join(output_path, organism_name, 'offtarget')

    # The stale middle layer of the chain.
    merged_path = os.path.join(offtarget_dir, f'{organism_name}_final_foldseek_results.json')
    with open(merged_path, 'w') as handle:
        json.dump({}, handle)

    offtargets.foldseek_human_parser(output_path, organism_name, mapping)
    merged = offtargets.merge_foldseek_data(output_path, organism_name)

    assert len(merged) == len(mapping), 'the stale empty merge was reused'
    assert sum(1 for value in merged.values() if value.get('target')) == 4


def test_the_final_table_is_not_reused_when_the_merge_changes(foldseek_results):
    output_path, organism_name, results_dir, mapping = foldseek_results
    offtarget_dir = os.path.join(output_path, organism_name, 'offtarget')

    # The stale last layer: a table that says No hit for every gene.
    table_path = os.path.join(offtarget_dir, f'{organism_name}_final_foldseek_results.tsv')
    with open(table_path, 'w') as handle:
        handle.write('gene\tFS_organism_structure_query\n')
        for gene in GENES:
            handle.write(f'{gene}\tNo hit\n')
    files.write_manifest(table_path, {'mapped_genes': 0, 'genes_with_hit': 0})

    offtargets.foldseek_human_parser(output_path, organism_name, mapping)
    merged = offtargets.merge_foldseek_data(output_path, organism_name)
    table = offtargets.final_foldseek_structure_table(output_path, organism_name, merged)

    assert len(table) == len(GENES), 'one row per gene regardless of hits'
    hits = table['FS_human_structure_hit'].notna().sum()
    assert hits == 4, f'the stale No-hit table was reused ({hits} hits)'


def test_the_whole_chain_runs_clean_from_nothing(foldseek_results):
    """With no stale files at all, the three steps agree on the same numbers."""
    output_path, organism_name, results_dir, mapping = foldseek_results

    parsed = offtargets.foldseek_human_parser(output_path, organism_name, mapping)
    merged = offtargets.merge_foldseek_data(output_path, organism_name)
    table = offtargets.final_foldseek_structure_table(output_path, organism_name, merged)

    parsed_hits = sum(1 for v in parsed.values() if v.get('target_foldseek'))
    merged_hits = sum(1 for v in merged.values() if v.get('target'))
    table_hits = int(table['FS_human_structure_hit'].notna().sum())

    assert parsed_hits == merged_hits == table_hits == 4


def test_a_gene_without_results_is_recorded_rather_than_dropped(foldseek_results):
    output_path, organism_name, results_dir, mapping = foldseek_results

    parsed = offtargets.foldseek_human_parser(output_path, organism_name, mapping)

    assert GENES[4] in parsed, 'a gene whose search produced nothing must still have an entry'
    assert parsed[GENES[4]]['target_foldseek'] is None


# ------------------------------- why there is no hit

def test_the_three_reasons_for_no_hit_are_distinguishable(foldseek_results):
    """
    'No human homolog' is a finding about the protein and a good sign for a drug
    target; 'the search never ran' is a failure of the run. They used to produce
    the same record with thirteen None fields.
    """
    output_path, organism_name, results_dir, mapping = foldseek_results

    # One gene searched with hits, one never searched, one with an unreadable file.
    unreadable = os.path.join(results_dir, 'AF_X_output_foldseek', 'AF_X_foldseek_results.tsv')
    os.makedirs(os.path.dirname(unreadable), exist_ok=True)
    open(unreadable, 'w').close()

    mapping = dict(mapping)
    mapping[GENES[4]] = []
    mapping['GENE_UNREADABLE'] = [unreadable]

    parsed = offtargets.foldseek_human_parser(output_path, organism_name, mapping)

    assert parsed[GENES[0]]['foldseek_status'] == 'hit'
    assert parsed[GENES[4]]['foldseek_status'] == 'not_searched'
    assert parsed['GENE_UNREADABLE']['foldseek_status'] == 'result_unreadable'


def test_a_search_with_no_hits_is_recorded_as_such(organism):
    output_path, organism_name = organism
    results_dir = os.path.join(output_path, organism_name, 'offtarget', 'foldseek_results')
    gene_dir = os.path.join(results_dir, 'AF_empty_output_foldseek')
    os.makedirs(gene_dir, exist_ok=True)
    path = os.path.join(gene_dir, 'AF_empty_foldseek_results.tsv')
    with open(path, 'w') as handle:
        # A header with no rows: Foldseek ran and found nothing.
        handle.write('query\ttarget\talnlen\tqcov\ttcov\tlddt\tqtmscore\tttmscore\t'
                     'alntmscore\trmsd\tprob\tpident\tevalue\n')

    parsed = offtargets.foldseek_human_parser(output_path, organism_name, {GENES[0]: [path]})

    assert parsed[GENES[0]]['foldseek_status'] == 'no_hits'
    assert parsed[GENES[0]]['target_foldseek'] is None


def test_the_reason_reaches_the_results_table(foldseek_results):
    output_path, organism_name, results_dir, mapping = foldseek_results

    offtargets.foldseek_human_parser(output_path, organism_name, mapping)
    merged = offtargets.merge_foldseek_data(output_path, organism_name)
    table = offtargets.final_foldseek_structure_table(output_path, organism_name, merged)

    assert 'FS_status' in table.columns
    statuses = dict(zip(table['gene'], table['FS_status']))
    assert statuses[GENES[0]] == 'hit'
    assert statuses[GENES[4]] == 'not_searched'


def test_a_gene_with_no_reference_structure_is_labelled_as_such(foldseek_results):
    """
    A gene the structures module never resolved is a third situation again: it
    was not searched because there was nothing to search with.
    """
    output_path, organism_name, results_dir, mapping = foldseek_results
    searched = {gene: paths for gene, paths in mapping.items() if paths}

    offtargets.foldseek_human_parser(output_path, organism_name, searched)
    merged = offtargets.merge_foldseek_data(output_path, organism_name)
    table = offtargets.final_foldseek_structure_table(output_path, organism_name, merged)

    statuses = dict(zip(table['gene'], table['FS_status']))
    assert statuses[GENES[4]] == 'no_reference_structure'
