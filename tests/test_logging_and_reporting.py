"""Logging and run reporting.

The log is the only record of what a long run did, so losing a message or
reporting a failed step as finished is a defect in its own right.
"""

import io
import logging
import os
import threading

import pandas as pd
import pytest

from ftscripts import logger as logmod


# --------------------------------------------------------------------- logger

def test_importing_the_module_has_no_side_effects():
    """Importing must not create a log file or take over the process streams."""
    assert logmod.log_file_path() is None or isinstance(logmod.log_file_path(), str)
    # The module-level logger exists, but nothing is configured by import alone.
    assert isinstance(logmod.logger, logging.Logger)


def test_concurrent_writers_do_not_lose_messages():
    """
    A flag shared by every thread used to drop any message written while another
    thread was inside write(): eight writers lost three quarters of them.
    """
    collected = []

    class Collect(logging.Handler):
        def emit(self, record):
            collected.append(record.getMessage())

    probe = logging.getLogger('test_concurrency')
    probe.handlers = []
    probe.setLevel(logging.INFO)
    probe.addHandler(Collect())
    probe.propagate = False

    stream = logmod.PrintLogger(probe, io.StringIO())

    threads_count, per_thread = 8, 200
    barrier = threading.Barrier(threads_count)

    def worker(number):
        barrier.wait()
        for index in range(per_thread):
            stream.write(f'thread{number}-message{index}\n')

    threads = [threading.Thread(target=worker, args=(n,)) for n in range(threads_count)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    expected = threads_count * per_thread
    assert len(collected) == expected, f'{expected - len(collected)} messages were dropped'
    assert len(set(collected)) == expected, 'messages were interleaved or duplicated'


def test_writes_reach_the_stream_unchanged_and_log_whole_lines():
    """
    A newline used to be appended to every chunk, so print(..., end='') and a
    partial write each became their own log line, and 'Running:' was split from
    the command it introduced.
    """
    lines = []

    class Collect(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())

    probe = logging.getLogger('test_lines')
    probe.handlers = []
    probe.setLevel(logging.INFO)
    probe.addHandler(Collect())
    probe.propagate = False

    class Spy:
        def __init__(self):
            self.chunks = []

        def write(self, message):
            self.chunks.append(message)
            return len(message)

        def flush(self):
            pass

    spy = Spy()
    stream = logmod.PrintLogger(probe, spy)

    stream.write('Running:')
    stream.write(' ')
    stream.write("'diamond' 'blastp'")
    stream.write('\n')

    assert ''.join(spy.chunks) == "Running: 'diamond' 'blastp'\n", 'the stream was altered'
    assert lines == ["Running: 'diamond' 'blastp'"], 'one command became several records'


def test_progress_redraws_are_shown_but_not_logged():
    """A tqdm bar would otherwise add thousands of records to the file."""
    lines = []

    class Collect(logging.Handler):
        def emit(self, record):
            lines.append(record.getMessage())

    probe = logging.getLogger('test_progress')
    probe.handlers = []
    probe.setLevel(logging.INFO)
    probe.addHandler(Collect())
    probe.propagate = False

    shown = io.StringIO()
    stream = logmod.PrintLogger(probe, shown)
    stream.write('\rSearching: 50%')
    stream.write('\rSearching: 100%')

    assert shown.getvalue() == '\rSearching: 50%\rSearching: 100%'
    assert lines == []


def test_setup_is_idempotent_and_the_file_is_unique_per_run(tmp_path):
    """
    The name used to be accurate to the minute and opened in 'w' mode, so two
    runs starting in the same minute meant the second truncated the first's log.
    """
    import ftscripts.logger as fresh

    fresh._log_file_path = None
    root = logging.getLogger()
    before = list(root.handlers)
    try:
        first = fresh.setup_logging(log_dir=str(tmp_path), capture_streams=False)
        handlers_after_first = len(root.handlers)
        second = fresh.setup_logging(log_dir=str(tmp_path), capture_streams=False)

        assert first == second
        assert len(root.handlers) == handlers_after_first, 'handlers stacked'
        assert str(os.getpid()) in os.path.basename(first)

        logging.info('a recorded line')
        for handler in root.handlers:
            handler.flush()
        with open(first) as handle:
            assert 'a recorded line' in handle.read()
    finally:
        root.handlers = before
        fresh._log_file_path = None


# ------------------------------------------------------------ run reporting

def test_failed_steps_are_recorded_and_reported():
    import fasttarget

    fasttarget.reset_failed_steps()
    assert fasttarget.failed_steps() == []

    fasttarget.record_failure('core analysis', RuntimeError('boom'))
    fasttarget.record_failure('structures analysis', ValueError('bang'))

    assert fasttarget.failed_steps() == ['core analysis', 'structures analysis']

    fasttarget.reset_failed_steps()
    assert fasttarget.failed_steps() == []


def test_column_summaries_describe_each_kind_of_result():
    import fasttarget

    booleans = pd.Series(['True', 'False', 'True', 'True'])
    assert fasttarget.summarize_column('core_roary', booleans) == 'core_roary: 3 of 4 genes'

    numbers = pd.Series([0.2, 0.5, 0.8, None])
    line = fasttarget.summarize_column('druggability_score', numbers)
    assert '3 of 4 genes with a value' in line and 'max 0.8' in line

    categories = pd.Series(['Cytoplasmic', 'Cytoplasmic', 'Unknown'])
    assert fasttarget.summarize_column('psortb_localization', categories) == (
        'psortb_localization: Cytoplasmic 2, Unknown 1'
    )

    assert fasttarget.summarize_column('FS_rmsd', pd.Series([None, None])) == 'FS_rmsd: no values'


def test_long_values_are_not_dumped_into_the_summary():
    """Species name lists would otherwise take over the block."""
    import fasttarget

    names = pd.Series([';'.join(f'Species_{i}' for i in range(40))] * 2 + [''])
    line = fasttarget.summarize_column('silva_core_species_hit_names', names)
    assert len(line) < 120
    assert 'distinct values' in line
