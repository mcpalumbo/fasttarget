"""The metabolism module: step reporting, chokepoint parsing and mapping.

Two themes recur here. A step must not report success it did not have, and an
input that is missing must not be replaced by the copy a previous run left in
the output directory.
"""

import os

import pytest

from conftest import GENES
from ftscripts import pathways


CHOKEPOINT_REPORT = """\
Reactions that are chokepoints on the producing side (3):
RXN-ONLY-PROD  produces something
RXN-BOTH-1  produces and consumes
RXN-BOTH-2  produces and consumes

Reactions that are chokepoints on the consuming side (3):
RXN-ONLY-CONS  consumes something
RXN-BOTH-1  produces and consumes
RXN-BOTH-2  produces and consumes

Reactions that are not chokepoints, but are candidates
RXN-NEITHER  nothing special
"""

EMPTY_SECTIONS_REPORT = """\
Reactions that are chokepoints on the producing side (0):

Reactions that are chokepoints on the consuming side (0):

Reactions that are not chokepoints, but are candidates
"""


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as handle:
        handle.write(text)
    return path


# ------------------------------------------------------------- step reporting

def test_a_step_that_succeeds_reports_finished(capsys):
    with pathways.metabolism_step(1, 'Doing something'):
        pass

    output = capsys.readouterr().out
    assert '1. Doing something' in output
    assert '1. Finished' in output
    assert 'FAILED' not in output


def test_a_step_that_fails_says_so_and_does_not_say_finished(capsys):
    """
    Every step used to print 'Finished' unconditionally, so a log read as a list
    of completed steps even where one had failed.
    """
    with pytest.raises(RuntimeError):
        with pathways.metabolism_step(3, 'Generating network SIF file'):
            raise RuntimeError('the tool exited 1')

    output = capsys.readouterr().out
    assert '3. FAILED: Generating network SIF file' in output
    assert '3. Finished' not in output


def test_the_original_exception_reaches_the_caller():
    """The failure has to be recordable by the caller, not swallowed here."""
    with pytest.raises(ValueError, match='specific cause'):
        with pathways.metabolism_step(2, 'Something'):
            raise ValueError('specific cause')


# ------------------------------------------------------- chokepoint parsing

def test_chokepoints_are_split_into_three_disjoint_sets(organism):
    output_path, organism_name = organism
    report = _write(os.path.join(str(output_path), 'chokepoints.txt'), CHOKEPOINT_REPORT)

    producing, consuming, both = pathways.process_chokepoint_ptools(
        output_path, organism_name, report
    )

    assert producing == {'RXN-ONLY-PROD'}
    assert consuming == {'RXN-ONLY-CONS'}
    assert both == {'RXN-BOTH-1', 'RXN-BOTH-2'}
    # Disjoint, which is what lets the mapping use an if/elif chain.
    assert not (producing & consuming) and not (producing & both) and not (consuming & both)


def test_an_empty_section_yields_no_chokepoints(organism, capsys):
    """
    ''.split('\\n') is [''], so an empty section used to produce one chokepoint
    named '' and the printed count reported it as real.
    """
    output_path, organism_name = organism
    report = _write(os.path.join(str(output_path), 'chokepoints.txt'), EMPTY_SECTIONS_REPORT)

    producing, consuming, both = pathways.process_chokepoint_ptools(
        output_path, organism_name, report
    )

    assert producing == set()
    assert consuming == set()
    assert both == set()
    assert '' not in producing | consuming | both

    printed = capsys.readouterr().out
    assert 'Total producing chokepoints\n0' in printed


def test_a_missing_chokepoint_file_does_not_fall_back_on_a_previous_run(organism):
    """
    The error used to be logged and the copy from an earlier run read instead, so
    a wrong path in the config produced another organism's chokepoints.
    """
    output_path, organism_name = organism
    stale = os.path.join(output_path, organism_name, 'metabolism',
                         f'{organism_name}_chokepoints.txt')
    _write(stale, CHOKEPOINT_REPORT)

    with pytest.raises(FileNotFoundError) as error:
        pathways.process_chokepoint_ptools(
            output_path, organism_name, os.path.join(str(output_path), 'absent.txt')
        )
    assert 'absent.txt' in str(error.value)


def test_a_report_without_the_expected_headings_is_empty_not_wrong(organism):
    output_path, organism_name = organism
    report = _write(os.path.join(str(output_path), 'chokepoints.txt'), 'unexpected content\n')

    producing, consuming, both = pathways.process_chokepoint_ptools(
        output_path, organism_name, report
    )

    assert (producing, consuming, both) == (set(), set(), set())


# ----------------------------------------------------------------- smarttable

SMARTTABLE = 'ID\tReactions\n' + '\n'.join([
    f'{GENES[0]}\tRXN-A, RXN-B',
    f'{GENES[1]}\tRXN-C',
    f'{GENES[2]}\t',
]) + '\n'


def test_the_smarttable_maps_genes_to_their_reactions(organism):
    output_path, organism_name = organism
    table = _write(os.path.join(str(output_path), 'smarttable.tsv'), SMARTTABLE)

    mapping = pathways.process_genes_smarttable(output_path, organism_name, table)

    assert mapping[GENES[0]] == ['RXN-A', 'RXN-B']
    assert mapping[GENES[1]] == ['RXN-C']
    assert GENES[2] not in mapping, 'a gene with no reactions is dropped, not empty'


def test_a_missing_smarttable_does_not_fall_back_on_a_previous_run(organism):
    output_path, organism_name = organism
    stale = os.path.join(output_path, organism_name, 'metabolism',
                         f'{organism_name}_smarttable.txt')
    _write(stale, SMARTTABLE)

    with pytest.raises(FileNotFoundError) as error:
        pathways.process_genes_smarttable(
            output_path, organism_name, os.path.join(str(output_path), 'absent.tsv')
        )
    assert 'absent.tsv' in str(error.value)


# -------------------------------------------------------------------- mapping

def test_a_gene_takes_the_highest_value_of_its_reactions():
    locus_reactions = {GENES[0]: ['RXN-A', 'RXN-B'], GENES[1]: ['RXN-C']}
    values = {'RXN-A': 0.2, 'RXN-B': 0.9, 'RXN-C': 0.5}

    mapped = pathways.mapping_metabolism(locus_reactions, values)

    assert mapped[GENES[0]] == 0.9
    assert mapped[GENES[1]] == 0.5


def test_a_gene_whose_reactions_are_absent_gets_zero_not_a_missing_key():
    locus_reactions = {GENES[0]: ['RXN-UNKNOWN'], GENES[1]: []}

    mapped = pathways.mapping_metabolism(locus_reactions, {'RXN-A': 0.5})

    assert mapped == {GENES[0]: 0, GENES[1]: 0}


def test_chokepoint_mapping_assigns_each_gene_to_one_category():
    locus_reactions = {
        GENES[0]: ['RXN-ONLY-PROD'],
        GENES[1]: ['RXN-ONLY-CONS'],
        GENES[2]: ['RXN-BOTH-1'],
        GENES[3]: ['RXN-NEITHER'],
    }

    producing, consuming, both = pathways.mapping_chokepoints(
        locus_reactions, {'RXN-ONLY-PROD'}, {'RXN-ONLY-CONS'}, {'RXN-BOTH-1'}
    )

    assert producing == {GENES[0]: 'RXN-ONLY-PROD'}
    assert consuming == {GENES[1]: 'RXN-ONLY-CONS'}
    assert both == {GENES[2]: 'RXN-BOTH-1'}
    assert GENES[3] not in producing and GENES[3] not in consuming and GENES[3] not in both


def test_a_gene_with_several_chokepoint_reactions_reports_the_last_one():
    """
    The table holds one value per gene, so only one reaction is reported. This
    pins that behaviour rather than endorsing it: a gene catalysing two producing
    chokepoints shows whichever comes last in its reaction list.
    """
    locus_reactions = {GENES[0]: ['RXN-P1', 'RXN-P2']}

    producing, _, _ = pathways.mapping_chokepoints(
        locus_reactions, {'RXN-P1', 'RXN-P2'}, set(), set()
    )

    assert producing == {GENES[0]: 'RXN-P2'}


# ------------------------------------------------- the stale SIF guard

def test_a_sif_that_was_not_rewritten_is_refused(organism, monkeypatch):
    """
    generate_sif reported success whether or not it had written anything, and the
    centrality step then read whichever network.sif was in the directory.
    """
    output_path, organism_name = organism
    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    _write(os.path.join(metabolism_dir, 'network.sif'), 'A\tinteracts\tB\n')

    monkeypatch.setattr(pathways.programs, 'run_sbml_to_sif',
                        lambda *args, **kwargs: None)

    with pytest.raises(RuntimeError) as error:
        pathways.generate_sif(output_path, organism_name)
    assert 'earlier run' in str(error.value)


def test_a_sif_that_was_written_is_accepted(organism, monkeypatch):
    output_path, organism_name = organism
    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    sif_path = os.path.join(metabolism_dir, 'network.sif')

    def write_sif(sbml_file, ubiquitous_file, out_dir):
        with open(os.path.join(out_dir, 'network.sif'), 'w') as handle:
            handle.write('A\tinteracts\tB\n')

    monkeypatch.setattr(pathways.programs, 'run_sbml_to_sif', write_sif)

    pathways.generate_sif(output_path, organism_name)

    assert os.path.exists(sif_path)


def test_a_missing_sif_is_refused(organism, monkeypatch):
    output_path, organism_name = organism

    monkeypatch.setattr(pathways.programs, 'run_sbml_to_sif',
                        lambda *args, **kwargs: None)

    with pytest.raises(FileNotFoundError):
        pathways.generate_sif(output_path, organism_name)


# ------------------------------------------------- ubiquitous compounds

def test_a_curated_file_is_used_when_given(organism):
    output_path, organism_name = organism
    curated = _write(os.path.join(str(output_path), 'curated.txt'), 'WATER\nATP\n')

    pathways.ubiquitous_checker(output_path, organism_name, curated)

    destination = os.path.join(output_path, organism_name, 'metabolism',
                               f'{organism_name}_ubiquitous.txt')
    with open(destination) as handle:
        assert handle.read() == 'WATER\nATP\n'


def test_a_curated_file_that_does_not_exist_is_an_error(organism):
    output_path, organism_name = organism

    with pytest.raises(FileNotFoundError):
        pathways.ubiquitous_checker(output_path, organism_name,
                                    os.path.join(str(output_path), 'absent.txt'))


def test_the_generated_file_is_used_when_no_curated_one_is_given(organism, capsys):
    output_path, organism_name = organism
    metabolism_dir = os.path.join(output_path, organism_name, 'metabolism')
    _write(os.path.join(metabolism_dir, 'ubiquitous_compounds.txt'), 'WATER\n')

    pathways.ubiquitous_checker(output_path, organism_name, None)

    destination = os.path.join(metabolism_dir, f'{organism_name}_ubiquitous.txt')
    assert os.path.exists(destination)
    assert 'auto-generated' in capsys.readouterr().out


def test_neither_file_is_an_error(organism):
    output_path, organism_name = organism

    with pytest.raises(FileNotFoundError):
        pathways.ubiquitous_checker(output_path, organism_name, None)
