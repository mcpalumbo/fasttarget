"""Logging setup for FastTarget.

The pipeline reports progress with `print`, and the log file is expected to hold
everything that reached the terminal. That is done by replacing `sys.stdout` and
`sys.stderr` with objects that both write through to the real stream and send
whole lines to a logger.

Setting that up is an explicit call, `setup_logging()`, not an import side
effect: importing a module should not create a directory, open a log file or
take over the process streams. `setup_logging()` is idempotent, so a module that
is imported twice, or a Nextflow task that calls it again, does not stack
handlers or reopen the file.
"""

import logging
import os
import sys
import threading
from datetime import datetime

# Kept so `from ftscripts.logger import logger` goes on working. Configuring the
# root logger means logging.info() and friends work from every module.
logger = logging.getLogger()

_setup_lock = threading.Lock()
_log_file_path = None


class PrintLogger:
    """A stream that mirrors what is written to it into a logger.

    Writes reach the terminal exactly as they were made, and complete lines are
    logged one record each. Partial writes are held until their newline arrives,
    so `print(..., end='')` and `sys.stdout.write()` of a fragment no longer each
    become a separate line in the log.
    """

    def __init__(self, target_logger, stream, level=logging.INFO):
        """
        :param target_logger: Logger receiving the completed lines.
        :param stream: Real stream to write through to.
        :param level: Level the lines are logged at.
        """
        self.logger = target_logger
        self.stream = stream
        self.level = level
        self._pending = ''
        self._lock = threading.RLock()
        self._local = threading.local()

    def _passthrough(self, message):
        self.stream.write(message)
        self.stream.flush()
        return len(message)

    def write(self, message):
        if not isinstance(message, str):
            message = str(message)

        # Progress bars rewrite the current line. Showing them is useful, logging
        # every redraw is not: one tqdm bar would otherwise add thousands of
        # records to the file.
        if message.startswith('\r'):
            return self._passthrough(message)

        # A logging handler writing to this stream would re-enter write() and
        # recurse. The previous version guarded that with a flag shared by every
        # thread, which silently discarded any message written while another
        # thread was inside write(): with eight concurrent writers three quarters
        # of them were lost, from the terminal as well as the log. The guard is
        # now per thread, and a re-entrant write still reaches the terminal
        # instead of disappearing.
        if getattr(self._local, 'logging_now', False):
            return self._passthrough(message)

        with self._lock:
            written = self._passthrough(message)

            self._pending += message
            lines = self._pending.split('\n')
            self._pending = lines.pop()

            self._local.logging_now = True
            try:
                for line in lines:
                    if line.strip():
                        # Only the right side: leading spaces carry the centring of
                        # a banner and the indentation of a traceback, which the
                        # previous strip() flattened.
                        self.logger.log(self.level, line.rstrip())
            finally:
                self._local.logging_now = False

        return written

    def flush(self):
        """Log whatever is held back and flush the real stream."""
        with self._lock:
            pending = self._pending
            self._pending = ''
            if pending.strip() and not getattr(self._local, 'logging_now', False):
                self._local.logging_now = True
                try:
                    self.logger.log(self.level, pending.rstrip())
                finally:
                    self._local.logging_now = False
        self.stream.flush()

    def isatty(self):
        """Report the real stream's answer, which tqdm and others check."""
        return getattr(self.stream, 'isatty', lambda: False)()

    @property
    def encoding(self):
        return getattr(self.stream, 'encoding', 'utf-8')


def log_file_path():
    """
    Path of the log file of this run.

    :return: Path of the log file, or None before setup_logging() has run.
    """
    return _log_file_path


def setup_logging(log_dir=None, level=logging.INFO, capture_streams=True):
    """
    Configure logging to a file and the console, once per process.

    :param log_dir: Directory for the log file. Defaults to 'logs' next to the
        repository root.
    :param level: Level for the file handler.
    :param capture_streams: Whether to route stdout and stderr into the log.
    :return: Path of the log file.
    """

    global _log_file_path

    with _setup_lock:
        if _log_file_path is not None:
            return _log_file_path

        if log_dir is None:
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            log_dir = os.path.join(parent_dir, 'logs')
        os.makedirs(log_dir, exist_ok=True)

        # The name used to be accurate to the minute and the file was opened in
        # 'w' mode, so two runs starting in the same minute, which is normal for
        # Nextflow tasks, meant the second truncated the first one's log. Seconds
        # and the pid make the name unique to the run.
        stamp = datetime.now().strftime('%Y-%m-%d-%H-%M-%S')
        _log_file_path = os.path.join(log_dir, f'fasttarget_{stamp}_{os.getpid()}.log')

        logger.setLevel(level)

        file_handler = logging.FileHandler(_log_file_path, mode='w')
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        logger.addHandler(file_handler)

        # Warnings and errors are shown on the console as well. The handler holds
        # the real stream, so logging never writes back through the replacements
        # installed below.
        console_handler = logging.StreamHandler(sys.__stderr__)
        console_handler.setLevel(logging.WARNING)
        console_handler.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        logger.addHandler(console_handler)

        if capture_streams:
            # Each replacement keeps its own real stream, instead of folding
            # stderr into stdout as before, so redirecting one of them still
            # behaves as expected.
            sys.stdout = PrintLogger(logger, sys.__stdout__)
            sys.stderr = PrintLogger(logger, sys.__stderr__)

        return _log_file_path
