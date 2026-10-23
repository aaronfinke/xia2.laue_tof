from __future__ import annotations

import contextlib
import logging
import pathlib
import shutil
import time

import xia2.Driver.timing

xia2_logger = logging.getLogger(__name__)

# Every program runs in PROCESSING; the files worth keeping are copied to RESULTS.
PROCESSING = "processing"
RESULTS = "results"


def report_timing(fn):
    def wrap_fn(*args, **kwargs):
        start_time = time.time()
        result = fn(*args, **kwargs)
        xia2_logger.debug("\nTiming report:")
        xia2_logger.debug("\n".join(xia2.Driver.timing.report()))
        duration = time.time() - start_time
        # write out the time taken in a human readable way
        xia2_logger.info(
            "Processing took %s", time.strftime("%Hh %Mm %Ss", time.gmtime(duration))
        )
        return result

    return wrap_fn


class _ResultFiles:
    """The logs and data files of a run, logs kept in the order they were written."""

    def __init__(self) -> None:
        self.logs: list[pathlib.Path] = []
        self.data: list[pathlib.Path] = []

    def record_log(self, path: pathlib.Path) -> None:
        if path.is_file() and path not in self.logs:
            self.logs.append(path)

    def record_data(self, path: pathlib.Path) -> None:
        if path.is_file() and path not in self.data:
            self.data.append(path)

    def collect(self, directory: pathlib.Path) -> None:
        """Copy everything recorded into directory, numbering the logs."""
        directory.mkdir(parents=True, exist_ok=True)
        for i, log in enumerate(self.logs, start=1):
            shutil.copyfile(log, directory / f"{i:02d}_{log.name}")
        for path in self.data:
            shutil.copyfile(path, directory / path.name)


result_files = _ResultFiles()


@contextlib.contextmanager
def collect_results(directory: pathlib.Path):
    """Copy the run's results into directory on the way out, even after a failure."""
    try:
        yield
    finally:
        result_files.collect(directory)
