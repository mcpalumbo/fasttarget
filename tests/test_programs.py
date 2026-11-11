"""Command assembly and subprocess reporting, with nothing executed.

These functions build the command lines the pipeline runs and report what
happened. Both are worth pinning without containers: a wrong argument silently
changes a result, and a failure that is not reported at ERROR cannot be found in
the log afterwards.
"""

import logging
import os
import subprocess

import pytest

from ftscripts import programs


@pytest.fixture
def containers(monkeypatch):
    """Record container invocations instead of running them."""
    calls = []

    def fake_run_container(work_dir, bind_dir, image_name, command, env_vars=None,
                           volumes=None, container_engine='docker'):
        calls.append({
            'work_dir': work_dir, 'bind_dir': bind_dir, 'image': image_name,
            'command': command, 'env_vars': env_vars, 'engine': container_engine,
        })

    monkeypatch.setattr(programs, 'run_container', fake_run_container)
    return calls


# --------------------------------------------------------- run_bash_command

def test_a_successful_command_returns_its_output():
    assert 'hello' in programs.run_bash_command(['echo', 'hello'])


def test_a_successful_command_is_recorded_at_debug_not_info(caplog):
    """
    Echoing every command at INFO put 4744 DIAMOND invocations and their output
    into one 23 MB log, which buried the lines describing what the run was doing.
    """
    with caplog.at_level(logging.DEBUG):
        programs.run_bash_command(['echo', 'quiet please'])

    levels = {record.levelno for record in caplog.records if 'quiet please' in record.message}
    assert levels and levels <= {logging.DEBUG}


def test_nothing_is_recorded_at_info_for_a_successful_command(caplog):
    with caplog.at_level(logging.INFO):
        programs.run_bash_command(['echo', 'hello'])

    assert [r for r in caplog.records if r.levelno >= logging.INFO] == []


def test_a_failing_command_raises():
    with pytest.raises(subprocess.CalledProcessError):
        programs.run_bash_command(['false'])


def test_a_failure_is_reported_at_error_with_the_command_and_its_stderr(caplog):
    """
    'Command failed' and the traceback used to be logged at INFO, so grepping the
    log for ERROR missed the cause of the failure.
    """
    with caplog.at_level(logging.DEBUG):
        with pytest.raises(subprocess.CalledProcessError):
            programs.run_bash_command(['sh', '-c', 'echo the reason >&2; exit 3'])

    errors = [r for r in caplog.records if r.levelno == logging.ERROR]
    assert errors, 'the failure was not reported at ERROR'
    message = errors[0].message
    assert 'return code 3' in message
    assert 'the reason' in message, 'stderr is what explains the failure'
    assert 'sh' in message, 'the command has to be identifiable'


def test_a_marker_is_written_for_a_completed_search(tmp_path):
    """
    An empty search result is valid, so the step records that it finished and a
    resume can tell it apart from a crash that left a zero-byte file.
    """
    from ftscripts import files

    output = os.path.join(str(tmp_path), 'search.tsv')
    open(output, 'w').close()

    assert files.step_completed(output) is False
    files.mark_complete(output)
    assert files.step_completed(output) is True


# ------------------------------------------------------- engine selection

@pytest.mark.parametrize('engine,expected', [
    ('singularity', True), ('apptainer', True), ('Singularity', True), ('APPTAINER', True),
    ('docker', False), ('Docker', False),
])
def test_singularity_like_engines_are_recognised_case_insensitively(engine, expected):
    assert programs._is_singularity_like_engine(engine) is expected


def test_an_unknown_engine_is_refused(monkeypatch):
    monkeypatch.setattr(programs, 'run_docker_container', lambda *a, **k: 'docker')
    monkeypatch.setattr(programs, 'run_singularity_container', lambda *a, **k: 'singularity')

    with pytest.raises(ValueError, match='podman'):
        programs.run_container('/work', '/bind', 'image', ['cmd'], container_engine='podman')


@pytest.mark.parametrize('engine,expected', [
    ('docker', 'docker'), ('singularity', 'singularity'), ('apptainer', 'singularity'),
])
def test_the_engine_decides_which_runner_is_used(monkeypatch, engine, expected):
    monkeypatch.setattr(programs, 'run_docker_container', lambda *a, **k: 'docker')
    monkeypatch.setattr(programs, 'run_singularity_container', lambda *a, **k: 'singularity')

    assert programs.run_container('/work', '/bind', 'image', ['cmd'],
                                  container_engine=engine) == expected


# --------------------------------------------------------- command assembly

def test_roary_receives_the_configured_thresholds(containers, tmp_path):
    """
    These were hardcoded to 99 and 95, so core.min_core_freq and core.min_identity
    had no effect on the run.
    """
    gff_dir = os.path.join(str(tmp_path), 'gff')
    os.makedirs(gff_dir)
    for name in ('a', 'b'):
        with open(os.path.join(gff_dir, f'{name}.gff'), 'w') as handle:
            handle.write('##gff-version 3\n')

    programs.run_roary(str(tmp_path), gff_dir, os.path.join(str(tmp_path), 'results'),
                       core_threshold=90, identity=80, cpus=3)

    command = containers[0]['command']
    assert command[command.index('-cd') + 1] == '90'
    assert command[command.index('-i') + 1] == '80'
    assert command[command.index('-p') + 1] == '3'


def test_roary_without_input_files_is_an_error(tmp_path, containers):
    gff_dir = os.path.join(str(tmp_path), 'gff')
    os.makedirs(gff_dir)

    with pytest.raises(Exception, match='No .gff files'):
        programs.run_roary(str(tmp_path), gff_dir, os.path.join(str(tmp_path), 'out'))


def test_corecruncher_receives_the_configured_thresholds(containers, tmp_path):
    faa_dir = os.path.join(str(tmp_path), 'faa')
    os.makedirs(faa_dir)
    with open(os.path.join(faa_dir, 'ref.faa'), 'w') as handle:
        handle.write('>a\nMKA\n')

    programs.run_core_cruncher(str(tmp_path), 'ref.faa', core_threshold=90, identity=80)

    command = containers[0]['command']
    assert command[command.index('-freq') + 1] == '90'
    assert command[command.index('-score') + 1] == '80'
    assert command[command.index('-ref') + 1] == 'ref.faa'


def test_fpocket_is_given_a_path_relative_to_the_work_directory(containers, tmp_path):
    """Absolute paths do not resolve inside a Singularity container."""
    pdb = os.path.join(str(tmp_path), 'models', 'AF_P1.pdb')
    os.makedirs(os.path.dirname(pdb))
    with open(pdb, 'w') as handle:
        handle.write('ATOM\n')

    programs.run_fpocket(str(tmp_path), pdb)

    command = containers[0]['command']
    assert command[0] == 'fpocket'
    assert command[-1] == os.path.join('models', 'AF_P1.pdb')
    assert not os.path.isabs(command[-1])


def test_psortb_gets_the_requested_gram_type(containers, tmp_path):
    # The input and the output directory are distinct, as in the pipeline, where
    # the FASTA is under genome/ and the results under localization/.
    faa = os.path.join(str(tmp_path), 'genome', 'proteins.faa')
    output_dir = os.path.join(str(tmp_path), 'localization')
    os.makedirs(os.path.dirname(faa))
    os.makedirs(output_dir)
    with open(faa, 'w') as handle:
        handle.write('>a\nMKA\n')

    programs.run_psort(faa, 'p', output_dir)

    assert '-p' in containers[0]['command']


@pytest.mark.parametrize('gram', ['x', 'N', 'gram-negative', ''])
def test_an_invalid_gram_type_is_refused(tmp_path, gram):
    faa = os.path.join(str(tmp_path), 'proteins.faa')
    with open(faa, 'w') as handle:
        handle.write('>a\nMKA\n')

    with pytest.raises(ValueError):
        programs.run_psort(faa, gram, str(tmp_path))


def test_an_invalid_psortb_output_format_is_refused(tmp_path):
    faa = os.path.join(str(tmp_path), 'proteins.faa')
    with open(faa, 'w') as handle:
        handle.write('>a\nMKA\n')

    with pytest.raises(ValueError):
        programs.run_psort(faa, 'n', str(tmp_path), output_format='fancy')


# ------------------------------------------------- the interpreter for SNDG

def test_the_sbml_steps_run_in_the_interpreter_of_the_current_environment(monkeypatch, tmp_path):
    """
    'python3' resolves through PATH to the base conda install, whose networkx 3.x
    no longer has read_gpickle, so these steps failed even though the
    environment's own networkx supports it.
    """
    import sys

    recorded = []
    monkeypatch.setattr(programs, 'run_bash_command', lambda cmd: recorded.append(cmd))

    sbml = os.path.join(str(tmp_path), 'model.sbml')
    with open(sbml, 'w') as handle:
        handle.write('<sbml/>')
    # run_ubiquitous checks its output, so give it something to find.
    with open(os.path.join(str(tmp_path), 'ubiquitous_compounds.txt'), 'w') as handle:
        handle.write('WATER\n')

    programs.run_ubiquitous(sbml, str(tmp_path))

    assert recorded[0][0] == sys.executable
    assert recorded[0][0] != 'python3'
    assert recorded[0][1:3] == ['-m', 'SNDG.Network.SBMLProcessor']


def test_a_missing_sbml_file_is_refused_rather_than_logged(tmp_path):
    with pytest.raises(FileNotFoundError):
        programs.run_ubiquitous(os.path.join(str(tmp_path), 'absent.sbml'), str(tmp_path))


def test_a_missing_output_directory_is_refused(tmp_path):
    sbml = os.path.join(str(tmp_path), 'model.sbml')
    with open(sbml, 'w') as handle:
        handle.write('<sbml/>')

    with pytest.raises(FileNotFoundError):
        programs.run_ubiquitous(sbml, os.path.join(str(tmp_path), 'absent_dir'))


def test_the_sif_step_also_uses_the_current_interpreter(monkeypatch, tmp_path):
    import sys

    recorded = []
    monkeypatch.setattr(programs, 'run_bash_command', lambda cmd: recorded.append(cmd))

    sbml = os.path.join(str(tmp_path), 'model.sbml')
    ubiquitous = os.path.join(str(tmp_path), 'ubiquitous.txt')
    for path in (sbml, ubiquitous):
        with open(path, 'w') as handle:
            handle.write('x\n')

    programs.run_sbml_to_sif(sbml, ubiquitous, str(tmp_path))

    assert recorded[0][0] == sys.executable
    assert '-f' in recorded[0]
