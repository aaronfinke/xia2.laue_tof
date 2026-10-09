"""
xia2.laue_tof: A processing pipeline for neutron Laue time-of-flight data,
using tools from the DIALS package.

The data are expected to be a single crystal measured as one or more stationary
exposures at different orientations, each exposure a stack of TOF bins. Each
orientation is processed independently: (bin) -> import -> find_spots -> index
-> refine -> tof_integrate.

The input has to be data dials.import can read, which for a time-of-flight
experiment means it has already been binned; raw event data must be reduced
before it is given to this pipeline, which stops with an error otherwise.
    xia2.laue_tof image=run_binned.h5 unit_cell=x space_group=y
Give one image= per crystal orientation, e.g.:
    xia2.laue_tof image=setting1.h5 image=setting2.h5 unit_cell=x space_group=y
An incident (vanadium) run and an empty-instrument run are optional but
recommended; they are passed on to dials.tof_integrate:
    xia2.laue_tof image=run.h5 vanadium_run=vanadium.nxs empty_run=empty.nxs

Refer to the individual DIALS program documentation for more details on
time-of-flight processing in DIALS.
"""

from __future__ import annotations

import logging
import pathlib
import sys
import traceback

import iotbx.phil
from dials.util.options import ArgumentParser

import xia2.Driver.timing
import xia2.Handlers.Streams
from xia2.Applications.xia2_main import write_citations
from xia2.Handlers.Citations import Citations
from xia2.Handlers.Files import cleanup
from xia2.Modules.Laue_TOF.laue_tof import full_phil_str, run_xia2_laue_tof

phil_scope = iotbx.phil.parse(full_phil_str)

xia2_logger = logging.getLogger(__name__)


def _reject_removed_parameters(args) -> None:
    """
    Stop on a top-level d_min= or d_max=, which used to cut the spot list.

    They were removed because they cannot work on TOF data, and leaving it at
    that would be worse than keeping them: phil resolves an unambiguous suffix,
    so `d_min=1.8` now quietly becomes `integration.calculated_d_min=1.8` and
    changes which reflections are integrated instead of failing. An old command
    line therefore has to be refused rather than reinterpreted.
    """
    for arg in args:
        name = arg.split("=", 1)[0].strip()
        if name in ("d_min", "d_max"):
            sys.exit(
                f"{name} has been removed. It set a resolution cut for"
                " spotfinding, which DIALS cannot apply to time-of-flight data:"
                " the cut is a pixel mask, and building one needs a fixed s0"
                " that a polychromatic beam does not have, so dials.find_spots"
                " fails with 'PolychromaticBeam has no fixed s0'.\nIf you meant"
                " the resolution to integrate to, that is"
                " integration.calculated_d_min, which this would otherwise have"
                " silently become."
            )


def run(args=sys.argv[1:]):
    """
    Parse the command line input, setup logging and run the Laue-TOF pipeline.
    """
    Citations.cite("dials-integration")
    Citations.cite("xia2")
    parser = ArgumentParser(
        usage="xia2.laue_tof image=run_binned.h5 unit_cell=x space_group=y",
        read_experiments=False,
        read_reflections=False,
        phil=phil_scope,
        check_format=False,
        epilog=__doc__,
    )
    _reject_removed_parameters(args)
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    xia2.Handlers.Streams.setup_logging(
        logfile="xia2.laue_tof.log", debugfile="xia2.laue_tof.debug.log"
    )
    # remove the xia2 handler from the dials logger.
    dials_logger = logging.getLogger("dials")
    dials_logger.handlers.clear()

    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        xia2_logger.info("The following parameters have been modified:\n%s", diff_phil)

    cwd = pathlib.Path.cwd()
    try:
        with cleanup(cwd):
            run_xia2_laue_tof(cwd, params)
    except ValueError as e:
        # Unlike xia2.ssx, which exits 0 here, a run that could not do what it
        # was asked exits non-zero, so that a caller can tell
        xia2_logger.error(f"Error: {e}")
        sys.exit(1)
    except FileNotFoundError as e:
        xia2_logger.error(e)
        sys.exit(1)
    except Exception as e:
        with (cwd / "xia2-error.txt").open(mode="w") as fh:
            traceback.print_exc(file=fh)
        xia2_logger.error("Error: %s", str(e))
        xia2_logger.info(traceback.format_exc())
        xia2_logger.warning(
            "Please send the contents of xia2.laue_tof.log and xia2-error.txt to xia2.support@gmail.com"
        )
        sys.exit(1)

    write_citations(program="xia2.laue_tof")
