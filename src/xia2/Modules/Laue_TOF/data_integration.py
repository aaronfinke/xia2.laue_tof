from __future__ import annotations

import logging
import math
import os
import pathlib
import shutil
import subprocess
from dataclasses import dataclass, field

from dials.array_family import flex
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dxtbx.model import Crystal, ExperimentList
from dxtbx.serialize import load

from xia2.Driver.timing import record_step
from xia2.Handlers.Files import FileHandler
from xia2.Handlers.Streams import banner
from xia2.Modules.Laue_TOF.laue_tof import (
    NO_FALLBACK,
    ExportParams,
    FileInput,
    IndexingParams,
    IntegrationParams,
    LaueTOFSetup,
    SpotfindingParams,
)

xia2_logger = logging.getLogger(__name__)

# Not a dials.index method: the label used when an orientation is supplied,
# which routes dials.index into the known orientation indexer instead.
KNOWN_ORIENTATION = "known_orientation"

# Foreground/background mask methods, as named by dials.tof_integrate
MASKS = ("ellipse", "seed_skewness")

# Below this many strong spots, what the spot size filter discarded is worth
# querying: a few hundred spots is enough to index but leaves nothing spare.
FEW_SPOTS = 500


@dataclass
class IndexingResult:
    """The outcome of one attempt at indexing an exposure."""

    method: str
    n_indexed: int = 0
    n_strong: int = 0
    rmsd_px: float | None = None
    crystal: Crystal | None = None
    accepted: bool = False
    reason: str = ""

    @property
    def fraction_indexed(self) -> float:
        return self.n_indexed / self.n_strong if self.n_strong else 0.0

    def summary(self) -> str:
        if not self.n_strong:
            return f"{self.method}: {self.reason}"
        rmsd = f"{self.rmsd_px:.2f} px" if self.rmsd_px is not None else "n/a"
        return (
            f"{self.method}: {self.n_indexed}/{self.n_strong} indexed"
            f" ({self.fraction_indexed:.1%}), RMSD {rmsd}"
            + ("" if self.accepted else f" - rejected, {self.reason}")
        )


@dataclass
class ExposureResult:
    """The outcome of processing one exposure of the run."""

    name: str
    image: str
    n_strong: int = 0
    indexing: IndexingResult | None = None
    integrated: bool = False
    attempts: list[IndexingResult] = field(default_factory=list)

    @property
    def indexed(self) -> bool:
        return self.indexing is not None and self.indexing.accepted


@dataclass
class IntegrationChoices:
    """
    The method and mask one run was integrated with.

    Every exposure of a run has to be integrated the same way for the batches
    to be comparable when they are scaled together, which is now automatic:
    they are integrated in one call. The choices are still carried in one place
    because the two-axis decision - which profile method survived, and which
    foreground mask won - is made once and reported once.
    """

    method: str | None = None
    mask: str | None = None


@dataclass
class IntegrationPass:
    """
    What one call of dials.tof_integrate is asked to integrate.

    A run makes either one pass over the reflections that were found, or that
    pass followed by a second over everything the geometry predicts. The two
    differ only in these three settings, so they travel together rather than
    being read from the parameters at the point of use.
    """

    integration_type: str = "observed"
    d_min: float | None = None
    wavelength_range: tuple[float, float] | None = None


def _executable(name: str, override: pathlib.Path | None = None) -> str:
    """Find a program, either where the user said it is or on $PATH."""
    if override:
        return os.fspath(override)
    found = shutil.which(name)
    if not found:
        raise ValueError(
            f"Unable to find {name} on $PATH. Either add it to $PATH or give its"
            f" location in the xia2.laue_tof phil scope."
        )
    return found


def _run_program(
    command: list[str], working_directory: pathlib.Path, name: str
) -> subprocess.CompletedProcess:
    """Run one of the external programs, raising on anything but success."""
    xia2_logger.debug("Running %s", " ".join(command))
    with record_step(name):
        result = subprocess.run(
            command, cwd=working_directory, capture_output=True, encoding="utf-8"
        )
    if result.returncode:
        raise ValueError(
            f"{name} returned error status:\n{result.stderr or result.stdout}"
        )
    return result


def _record_log(tag: str, logfile: pathlib.Path) -> None:
    if logfile.is_file():
        FileHandler.record_log_file(tag, str(logfile))


# Neutron wavelength from time of flight: lambda[A] = H_OVER_M * t[s] / L[m]
H_OVER_M = 3956.034


def wavelength_range(experiment) -> tuple[float, float]:
    """
    The wavelength range of an experiment, in Angstrom, from its TOF range.

    Derived from the scan rather than taken from beam.get_wavelength_range(),
    because a format class may state a fixed band instead of the one the data
    holds: FormatESSNMX derives it (and the two agree to three decimals) while
    FormatMANDI returns a hard-coded 2.0-4.0 where the histogram gives
    1.92-4.08. Falls back to the beam if the scan has no time_of_flight.
    """
    scan = experiment.scan
    if scan is None or not scan.has_property("time_of_flight"):
        return tuple(experiment.beam.get_wavelength_range())  # type: ignore[return-value]
    tof = scan.get_property("time_of_flight")  # microseconds
    t_min, t_max = min(tof) * 1e-6, max(tof) * 1e-6
    to_source = experiment.beam.get_sample_to_source_distance() / 1000.0  # m
    distances = [panel.get_distance() / 1000.0 for panel in experiment.detector]
    longest = to_source + max(distances)
    shortest = to_source + min(distances)
    return H_OVER_M * t_min / longest, H_OVER_M * t_max / shortest


def geometric_d_min(experiment) -> tuple[float, float]:
    """
    The best resolution the detector and the wavelength band can reach, and the
    largest scattering angle, in Angstrom and degrees.

    d = lambda / (2 sin theta), so the limit is the shortest wavelength at the
    largest scattering angle any pixel sees. It is a property of the instrument
    and the binning, not of the crystal or of what was observed, which makes it
    the sensible place to stop when integrating every predicted reflection
    rather than only the spots that were found.
    """
    lambda_min, _ = wavelength_range(experiment)
    to_source = experiment.beam.get_sample_to_source_distance() / 1000.0
    unit_s0 = experiment.beam.get_unit_s0()
    best_d_min = None
    widest = 0.0
    for panel in experiment.detector:
        nx, ny = panel.get_image_size()
        corners = ((0, 0), (nx, 0), (0, ny), (nx, ny), (nx / 2, ny / 2))
        two_theta = max(
            panel.get_two_theta_at_pixel(unit_s0, corner) for corner in corners
        )
        if two_theta <= 0:
            continue
        # This panel's own flight path sets its shortest wavelength
        panel_lambda_min = lambda_min
        if experiment.scan is not None and experiment.scan.has_property(
            "time_of_flight"
        ):
            t_min = min(experiment.scan.get_property("time_of_flight")) * 1e-6
            panel_lambda_min = (
                H_OVER_M * t_min / (to_source + panel.get_distance() / 1000.0)
            )
        d_min = panel_lambda_min / (2 * math.sin(two_theta / 2))
        best_d_min = d_min if best_d_min is None else min(best_d_min, d_min)
        widest = max(widest, two_theta)
    if best_d_min is None:
        raise ValueError("No panel of this detector sees a non-zero scattering angle")
    return best_d_min, math.degrees(widest)


def _format_class_name(image: str) -> str:
    """
    The dxtbx format class that will read a file, for error messages.

    Which class claims a file matters here: a reduced MANDI file and a reduced
    NMX file are both NXlauetof, and an instrument-specific class that claims
    the wrong one of them fails deep inside dxtbx, where the traceback says
    nothing about why.
    """
    try:
        from dxtbx.format.Registry import get_format_class_for_file

        format_class = get_format_class_for_file(image)
    except Exception:
        return "unknown"
    return format_class.__name__ if format_class else "none"


def apply_detector_calibration(experiments: pathlib.Path, detcal: pathlib.Path) -> None:
    """
    Replace the imported detector geometry with an ISAW/Mantid .DetCal one.

    The nominal MANDI geometry places individual panels several pixels out.
    Refinement does not correct that: with one detector group and many panels
    it absorbs the error into the crystal orientation instead, which leaves
    per-panel residuals that no mask or profile model can recover from.  On
    13378 the calibration takes indexing from 222/329 reflections at a mean
    residual of 3.28 px to 314/329 at 1.38 px, and agreement with Mantid's own
    peak positions from 28.5 % to 74.4 % within 3 px.

    Conventions: DetCal Base -> fast axis, Up -> slow axis, Center is the panel
    midpoint in cm, and DETNUM matches the dxtbx panel named "bank<DETNUM>".
    Base and Up are written to five decimals, so slow is re-orthogonalised
    against fast before use.
    """
    import numpy as np

    cal = {}
    for line in detcal.read_text().splitlines():
        f = line.split()
        if len(f) >= 17 and f[0] == "5":
            cal[int(f[1])] = (
                float(f[4]) * 10.0,
                float(f[5]) * 10.0,
                np.array([float(f[8]), float(f[9]), float(f[10])]) * 10.0,
                np.array([float(f[11]), float(f[12]), float(f[13])]),
                np.array([float(f[14]), float(f[15]), float(f[16])]),
            )
    if not cal:
        raise ValueError(f"No panel records (lines starting '5') in {detcal}")

    expts = load.experiment_list(experiments, check_format=False)
    detector = expts[0].detector
    applied, missing, mismatched = 0, [], []
    for panel in detector:
        digits = "".join(c for c in panel.get_name() if c.isdigit())
        if not digits or int(digits) not in cal:
            missing.append(panel.get_name())
            continue
        width, height, centre, base, up = cal[int(digits)]
        fast = base / np.linalg.norm(base)
        slow = up - np.dot(up, fast) * fast
        slow /= np.linalg.norm(slow)
        if (
            float(np.dot(fast, np.array(panel.get_fast_axis()))) < 0.9
            or float(np.dot(slow, np.array(panel.get_slow_axis()))) < 0.9
        ):
            mismatched.append(panel.get_name())
            continue
        origin = centre - 0.5 * width * fast - 0.5 * height * slow
        panel.set_frame(tuple(fast), tuple(slow), tuple(origin))
        applied += 1
    if mismatched:
        raise ValueError(
            f"{detcal.name} does not describe this detector: the axes of"
            f" {len(mismatched)} panel(s) disagree by more than 25 deg, starting"
            f" with {mismatched[0]}"
        )
    expts.as_file(experiments)
    xia2_logger.info(
        f"Applied detector calibration {detcal.name} to {applied} panel(s)"
        + (f", {len(missing)} not in the file" if missing else "")
    )


def run_import(
    working_directory: pathlib.Path,
    inputs: list[tuple[str, str]],
    file_input: FileInput,
) -> pathlib.Path:
    """
    Run dials.import on every exposure of the run at once, and check that the
    data really is TOF.

    All the exposures go into one experiment list, which is what lets the
    following steps see the run as a whole. Each becomes its own experiment and
    its own imageset, because they are separate stationary measurements.
    """
    experiments = working_directory / "imported.expt"
    command = [
        _executable("dials.import"),
        "output.experiments=imported.expt",
    ]
    for path_type, path in inputs:
        command.append(f"{path_type}={path}" if path_type != "image" else path)
    if file_input.import_phil:
        command.insert(1, os.fspath(file_input.import_phil))
    if file_input.mask:
        command.append(f"mask={os.fspath(file_input.mask)}")

    xia2_logger.notice(banner("Importing"))  # type: ignore
    image = inputs[0][1]
    try:
        _run_program(command, working_directory, "dials.import")
    except ValueError as e:
        raise ValueError(
            f"dials.import failed on {image}, which dxtbx reads with"
            f" {_format_class_name(image)}.\n{e}"
        ) from e
    _record_log(
        f"{working_directory.name} import", working_directory / "dials.import.log"
    )

    if file_input.detector_calibration:
        apply_detector_calibration(experiments, file_input.detector_calibration)

    expts = load.experiment_list(experiments, check_format=False)
    if not expts.all_tof():
        raise ValueError(
            f"{image} was not imported as time-of-flight data. It was read with"
            f" {_format_class_name(image)}; check that the file holds TOF bins,"
            " and that the expected format class is reading it."
        )
    if len(expts) != len(inputs):
        raise ValueError(
            f"{len(inputs)} input(s) were given but dials.import made"
            f" {len(expts)} experiment(s). Every exposure of a run has to import"
            " as exactly one experiment for the batches to line up with the"
            " inputs."
        )
    _check_uniform(expts, inputs)
    xia2_logger.info(
        f"Imported {len(expts)} exposure(s) with {_format_class_name(image)},"
        f" {len(expts[0].detector)} panels, {expts[0].scan.get_num_images()} TOF bins"
    )
    low, high = wavelength_range(expts[0])
    d_min, two_theta = geometric_d_min(expts[0])
    xia2_logger.info(
        f"Wavelengths {low:.2f}-{high:.2f} A over scattering angles up to"
        f" {two_theta:.0f} deg, so the geometry allows d_min = {d_min:.3f} A"
    )
    return experiments


def _check_uniform(experiments, inputs: list[tuple[str, str]]) -> None:
    """
    Require every exposure to share a detector shape and a TOF binning.

    dials.tof_integrate reads the panel size and the number of TOF bins from
    the first experiment and applies them to all of them
    (`image_size = experiments[0].detector[0].get_image_size()`), so a run whose
    exposures were reduced to different bin counts is integrated against the
    wrong frame axis rather than rejected. Two reductions of the same NMX run,
    one at 50 bins and one at 100, import together without complaint, which is
    how this was found.
    """
    shapes: dict[tuple[int, int], list[str]] = {}
    for expt, (_, path) in zip(experiments, inputs):
        key = (len(expt.detector), expt.scan.get_num_images())
        shapes.setdefault(key, []).append(path)
    if len(shapes) > 1:
        lines = [
            f"  {panels} panels, {bins} TOF bins: {', '.join(paths)}"
            for (panels, bins), paths in sorted(shapes.items())
        ]
        raise ValueError(
            "The exposures of a run do not share a detector shape and a TOF"
            " binning, so they cannot be integrated together:\n"
            + "\n".join(lines)
            + "\nReduce them all to the same number of bins, or process them as"
            " separate runs."
        )


def _report_size_filter(
    logfile: pathlib.Path, params: SpotfindingParams, n_strong: int
) -> None:
    """
    Say how much of the spot list the minimum spot size took.

    min_spot_size is the parameter most worth revisiting on a new instrument:
    the default suits NMX, where small split peaks are the problem, but on data
    whose spots span few TOF frames it can discard most of what was found. The
    counts are in the dials.find_spots log, which reports what it extracted
    before filtering.

    Advice is only worth giving when the spots that survive are few: a
    thresholding pass turns up tens of thousands of one and two pixel
    candidates, so on NMX 99% of them are filtered out and 2296 spots remain,
    which is not a problem. On MANDI it is 59% of 558, leaving 227.
    """
    if not logfile.is_file():
        return
    extracted = removed = 0
    for line in logfile.read_text().splitlines():
        if line.startswith("Extracted ") and line.endswith(" spots"):
            extracted = int(line.split()[1])
        elif line.startswith("Removed ") and " with size < " in line:
            removed = int(line.split()[1])
    if not extracted or not removed:
        return
    fraction = removed / extracted
    message = (
        f"{removed} of {extracted} spots found were smaller than"
        f" spotfinding.min_spot_size ({params.min_spot_size} pixels)"
    )
    if fraction >= 0.5 and params.min_spot_size > 3 and n_strong < FEW_SPOTS:
        # Nothing to suggest once the size filter is down to a few pixels: a
        # threshold algorithm like radial_profile finds tens of thousands of
        # one and two pixel candidates by design, and dropping those is the
        # filter working rather than a loss.
        xia2_logger.warning(
            f"{message}. That is most of them: try a smaller"
            " spotfinding.min_spot_size, and"
            " spotfinding.threshold_algorithm=radial_profile, if indexing"
            " struggles for want of spots."
        )
    else:
        xia2_logger.info(message)


def find_spots(
    working_directory: pathlib.Path, params: SpotfindingParams, n_exposures: int = 1
) -> list[int]:
    """
    Run dials.find_spots over the whole run, returning the spots per exposure.

    One call covers every exposure, and the strong spots come back in a single
    table whose `id` column says which exposure each spot belongs to.

    There is deliberately no resolution cut here. `spotfinder.filter.d_min` and
    `d_max` cannot work on TOF data: both route through
    `dials/util/masking.py::_apply_resolution_mask`, which builds a
    `ResolutionMaskGenerator` from the beam, and a `PolychromaticBeam` has no
    fixed `s0` to build it from - `dials.find_spots` dies with
    `PolychromaticBeam has no fixed s0` at `beam.h(741)`. A pixel mask is the
    wrong mechanism anyway, since on a Laue-TOF image the resolution of a pixel
    depends on the wavelength of the neutron that arrived there. A cut would
    have to be made in reflection space, after the wavelength is known.
    """
    command = [
        _executable("dials.find_spots"),
        "imported.expt",
        "output.reflections=strong.refl",
        f"spotfinder.filter.min_spot_size={params.min_spot_size}",
        f"spotfinder.threshold.algorithm={params.threshold_algorithm}",
        f"spotfinder.mp.nproc={params.nproc}",
    ]
    if params.phil:
        command.insert(1, os.fspath(params.phil))
    if params.max_spot_size:
        command.append(f"spotfinder.filter.max_spot_size={params.max_spot_size}")

    xia2_logger.notice(banner("Spotfinding"))  # type: ignore
    _run_program(command, working_directory, "dials.find_spots")
    _record_log(
        f"{working_directory.name} find_spots",
        working_directory / "dials.find_spots.log",
    )

    strong = flex.reflection_table.from_file(working_directory / "strong.refl")
    counts = [(strong["id"] == i).count(True) for i in range(n_exposures)]
    total = strong.size()
    if n_exposures > 1:
        xia2_logger.info(
            f"Found {total} strong spots over {n_exposures} exposures: "
            + ", ".join(f"{n}" for n in counts)
        )
    else:
        xia2_logger.info(f"Found {total} strong spots")
    _report_size_filter(working_directory / "dials.find_spots.log", params, total)

    # The limit is per exposure, not over the run. What it guards against is an
    # assertion in the DIALS TOF reflection predictor, and refinement treats
    # each exposure as its own job - dials.refine says "independently-refined
    # groups of experiments have been recombined" - so one exposure's spot list
    # is what reaches the predictor, however many exposures there are.
    over = [i for i, n in enumerate(counts) if n > params.max_strong]
    if over:
        worst = max(over, key=lambda i: counts[i])
        raise ValueError(
            f"Exposure {worst + 1} has {counts[worst]} strong spots, more than"
            f" spotfinding.max_strong ({params.max_strong})"
            + (f"; {len(over)} exposures are over the limit" if len(over) > 1 else "")
            + ". Large spot lists trip an assertion in the DIALS TOF reflection"
            " predictor during refinement; raise spotfinding.min_spot_size, or"
            " raise the limit if you are sure."
        )
    return counts


def _primitive_setting_a_matrix(crystal: Crystal) -> tuple[float, ...]:
    """
    The A matrix of a crystal, in the primitive setting. Needed currently
    due to a DIALS bug.

    dials.index expects known_symmetry.A_matrix in the primitive setting and
    contains the conversion, but IndexerKnownOrientation assigns it to the loop
    variable and discards it (DIALS 3.30), while the post-indexing conversion
    back to the centred setting still runs - so a centred cell comes back wrong.
    Converting here supplies the setting the indexer needs, leaving the skipped
    conversion nothing to do. Note this is tied to the defect: the crystal is
    built with known_symmetry.space_group, so a fixed DIALS would convert this
    matrix a second time. _assess_indexing's cell check is what would catch it.
    """
    cb_op = crystal.get_space_group().info().change_of_basis_op_to_primitive_setting()
    return crystal.change_basis(cb_op).get_A()


def _assess_indexing(
    method: str,
    params: IndexingParams,
    n_strong: int,
    crystal: Crystal | None,
    indexed: flex.reflection_table,
) -> IndexingResult:
    """
    Decide whether one exposure's indexing solution is good enough to keep.

    Indexing that 'worked' is not enough: a small, clean spot list indexes
    almost completely, and an alternative indexing can index plenty of
    reflections on the wrong cell. So check the cell, the number indexed and
    the positional residuals.

    Takes the crystal and the exposure's own indexed reflections rather than a
    pair of file names, because one dials.index call now covers the whole run
    and each exposure has to be judged on its own rows.
    """
    result = IndexingResult(method=method, n_strong=n_strong)
    result.n_indexed = indexed.size()

    if "xyzcal.px" in indexed and indexed.size():
        # The median residual, not the mean: a solution that is right about most
        # of the spot list can carry a tail of badly fitted reflections that the
        # mean is dominated by, and rejecting it loses a whole orientation. On
        # MANDI 13380 the mean says 9.38 px and the median 1.92 px, where
        # dials.index, which reports its RMSD after refinement has discarded
        # outliers, says 1.45 / 1.12 px. Comparable with the dials.index log to
        # within the difference between a median and an outlier-rejected mean.
        dx, dy, _ = (indexed["xyzobs.px.value"] - indexed["xyzcal.px"]).parts()
        squared = flex.pow2(dx) + flex.pow2(dy)
        result.rmsd_px = float(sorted(squared)[squared.size() // 2] ** 0.5)

    if crystal is None:
        result.reason = "no crystal model"
        return result
    result.crystal = crystal

    required = max(params.min_indexed, int(params.min_indexed_fraction * n_strong))
    if result.n_indexed < required:
        result.reason = (
            f"only {result.n_indexed} of {n_strong} strong spots indexed"
            f" ({result.fraction_indexed:.1%}), fewer than the {required}"
            f" required by indexing.min_indexed ({params.min_indexed}) and"
            f" indexing.min_indexed_fraction ({params.min_indexed_fraction:.0%})"
        )
        return result

    if (
        params.max_rmsd_px
        and result.rmsd_px is not None
        and result.rmsd_px > params.max_rmsd_px
    ):
        result.reason = (
            f"RMSD {result.rmsd_px:.2f} px is worse than"
            f" indexing.max_rmsd_px ({params.max_rmsd_px})"
        )
        return result

    if params.unit_cell and not crystal.get_unit_cell().is_similar_to(
        params.unit_cell, relative_length_tolerance=0.1, absolute_angle_tolerance=5.0
    ):
        result.reason = (
            f"refined cell {crystal.get_unit_cell()} does not match the"
            f" target cell {params.unit_cell}"
        )
        return result

    result.accepted = True
    return result


def _write_indexing_subset(
    working_directory: pathlib.Path, exposures: list[int], tag: str
) -> tuple[str, str]:
    """
    Write the imported experiments and strong spots of a subset of exposures.

    A rung of the ladder only has to re-index the exposures the earlier rungs
    did not solve, so each rung gets its own pair of input files holding just
    those. Experiment ids are renumbered from zero, which is what dials.index
    expects, and the mapping back to the original exposure is the position in
    `exposures`.
    """
    expts = load.experiment_list(
        working_directory / "imported.expt", check_format=False
    )
    strong = flex.reflection_table.from_file(working_directory / "strong.refl")

    subset = ExperimentList([expts[i] for i in exposures])
    tables = []
    for new_id, old_id in enumerate(exposures):
        rows = strong.select(strong["id"] == old_id)
        rows["id"] = flex.int(rows.size(), new_id)
        if "imageset_id" in rows:
            rows["imageset_id"] = flex.int(rows.size(), new_id)
        tables.append(rows)
    table = flex.reflection_table.concat(tables) if len(tables) > 1 else tables[0]

    experiments_file = f"to_index_{tag}.expt"
    reflections_file = f"to_index_{tag}.refl"
    subset.as_file(working_directory / experiments_file)  # type: ignore[attr-defined]
    table.as_file(working_directory / reflections_file)
    return experiments_file, reflections_file


def _indexing_base_command(params: IndexingParams) -> list[str]:
    """The dials.index arguments every rung of the ladder shares."""
    command = [_executable("dials.index")]
    if params.phil:
        command.append(os.fspath(params.phil))
    if params.space_group:
        command.append(f"indexing.known_symmetry.space_group={params.space_group}")
    if params.unit_cell:
        cell = ",".join(f"{p}" for p in params.unit_cell.parameters())
        command.append(f"indexing.known_symmetry.unit_cell={cell}")
    if params.max_cell:
        command.append(f"indexing.max_cell={params.max_cell}")
    if params.d_min_start:
        command.append(f"indexing.refinement_protocol.d_min_start={params.d_min_start}")
    command += [
        # The exposures of a run are separate stationary measurements of one
        # crystal, each at its own setting, so there is no single orientation to
        # fit and they must be indexed independently. DIALS cannot infer this:
        # a TOF experiment carries a goniometer and a scan over its bins, so
        # is_still() is False and joint_indexing=Auto takes the all-rotation
        # branch, which refuses to choose and aborts the run.
        "indexing.joint_indexing=False",
        f"refinement.reflections.outlier.algorithm={params.outlier_algorithm}",
        f"refinement.reflections.outlier.tukey.iqr_multiplier={params.outlier_iqr_multiplier}",
        f"refinement.reflections.outlier.separate_panels={params.outlier_separate_panels}",
    ]
    return command


def _assess_rung(
    working_directory: pathlib.Path,
    tag: str,
    method: str,
    params: IndexingParams,
    exposures: list[int],
    n_strong: list[int],
) -> tuple[dict[int, IndexingResult], dict[int, int]]:
    """
    Judge what one dials.index call returned, per exposure.

    The output holds only the exposures that produced a crystal, in the order
    their worker processes finished - dials.index extends the list inside an
    `as_completed` loop - so the position of an experiment says nothing about
    which exposure it came from. Its reflections' `imageset_id` does: DIALS
    resets that to the index within this call, which maps back through
    `exposures`.

    Returns the assessment per exposure, and where each exposure's experiment
    sits in this rung's output, which is what lifts it back out later.
    """
    experiments_file = working_directory / f"indexed_{tag}.expt"
    reflections_file = working_directory / f"indexed_{tag}.refl"
    results: dict[int, IndexingResult] = {}
    positions: dict[int, int] = {}
    if not (experiments_file.is_file() and reflections_file.is_file()):
        return results, positions

    expts = load.experiment_list(experiments_file, check_format=False)
    refl = flex.reflection_table.from_file(reflections_file)
    # A rung that solved nothing still writes its files, and the reflection
    # table is then not merely empty but column-less, so asking it for a flag
    # raises rather than returning nothing.
    if not refl.size() or "flags" not in refl:
        return results, positions

    indexed = refl.select(refl.get_flags(refl.flags.indexed))
    for position, expt in enumerate(expts):
        rows = indexed.select(indexed["id"] == position)
        if not rows.size():
            continue
        subset_ids = set(rows["imageset_id"])
        if len(subset_ids) != 1:
            raise ValueError(
                f"Experiment {position} of {experiments_file.name} mixes"
                f" imageset ids {sorted(subset_ids)}, so it cannot be traced"
                " back to one exposure"
            )
        exposure = exposures[subset_ids.pop()]
        positions[exposure] = position
        results[exposure] = _assess_indexing(
            method, params, n_strong[exposure], expt.crystal, rows
        )
    return results, positions


def _merge_indexed(
    working_directory: pathlib.Path,
    chosen: dict[int, tuple[str, int]],
) -> None:
    """
    Gather the accepted solutions into one indexed.expt / indexed.refl.

    The ladder leaves each rung's output on disk, so an accepted exposure has
    to be lifted out of whichever rung solved it. Exposures are written in
    input order, so the batches that reach the MTZ follow the order the
    exposures were given in rather than the order they happened to index in.

    Only the indexed reflections are carried over. The unindexed ones share a
    single id of -1, which says nothing about which exposure they came from
    once several exposures are in one table, and nothing downstream reads them.
    """
    experiments = ExperimentList()
    tables = []
    for exposure in sorted(chosen):
        tag, position = chosen[exposure]
        expts = load.experiment_list(
            working_directory / f"indexed_{tag}.expt", check_format=False
        )
        refl = flex.reflection_table.from_file(
            working_directory / f"indexed_{tag}.refl"
        )
        rows = refl.select(refl["id"] == position)
        new_id = len(experiments)
        rows["id"] = flex.int(rows.size(), new_id)
        rows["imageset_id"] = flex.int(rows.size(), new_id)
        experiments.append(expts[position])
        tables.append(rows)

    table = flex.reflection_table.concat(tables) if len(tables) > 1 else tables[0]
    experiments.as_file(working_directory / "indexed.expt")  # type: ignore[attr-defined]
    table.as_file(working_directory / "indexed.refl")


def index(
    working_directory: pathlib.Path,
    params: IndexingParams,
    n_strong: list[int],
) -> list[list[IndexingResult]]:
    """
    Index every exposure of the run, walking the method ladder over the set.

    One dials.index call covers all the exposures still wanting a solution, with
    joint_indexing=False so that each is indexed on its own. Near the limit of
    what a detector layout can support indexing succeeds stochastically rather
    than reproducibly, so a rung that leaves some exposures unsolved is followed
    by the next method on just those, and finally by a seeded attempt using the
    crystal of an exposure that did solve - which is what rescues layouts that
    never index de novo.

    Returns the attempts made for each exposure, the last of which is its
    accepted solution if it has one.
    """
    attempts: list[list[IndexingResult]] = [[] for _ in n_strong]
    pending = []
    for exposure, count in enumerate(n_strong):
        if count < params.min_spots:
            attempts[exposure].append(
                IndexingResult(
                    method="none",
                    n_strong=count,
                    reason=f"only {count} strong spots, fewer than indexing.min_spots",
                )
            )
        else:
            pending.append(exposure)

    if not pending:
        return attempts

    base_command = _indexing_base_command(params)
    xia2_logger.notice(banner("Indexing"))  # type: ignore

    # Where every accepted solution lives, keyed by exposure and method: the
    # rung's tag and the position of that exposure's experiment in the rung's
    # output. Keyed by method as well as exposure because an exposure can be
    # solved more than once, and the best of those is chosen at the end.
    locations: dict[tuple[int, str], tuple[str, int]] = {}
    # The best accepted result per exposure so far, to pick a seed crystal from.
    accepted: dict[int, IndexingResult] = {}

    rungs: list[tuple[str, str]] = [(method, method) for method in params.methods]
    if params.seed_from_first:
        rungs.append((KNOWN_ORIENTATION, KNOWN_ORIENTATION))

    for method, tag in rungs:
        if not pending:
            break
        seed_crystal: Crystal | None = None
        if method == KNOWN_ORIENTATION:
            solved = [
                result for result in accepted.values() if result.crystal is not None
            ]
            if not solved:
                # Nothing indexed, so there is no model to seed with.
                break
            best = min(solved, key=lambda r: (r.rmsd_px is None, r.rmsd_px))
            seed_crystal = best.crystal
            xia2_logger.info(
                f"Seeding the {len(pending)} unindexed exposure(s) with the"
                " crystal of the best solution so far"
            )

        experiments_file, reflections_file = _write_indexing_subset(
            working_directory, pending, tag
        )
        command = base_command + [
            experiments_file,
            reflections_file,
            f"output.log=dials.index.{tag}.log",
            f"output.experiments=indexed_{tag}.expt",
            f"output.reflections=indexed_{tag}.refl",
        ]
        if seed_crystal is not None:
            # Supplying an A matrix routes dials.index into the known
            # orientation indexer, which ignores indexing.method.
            a_matrix = ",".join(
                f"{v}" for v in _primitive_setting_a_matrix(seed_crystal)
            )
            command.append(f"indexing.known_symmetry.A_matrix={a_matrix}")
        else:
            command.append(f"indexing.method={method}")

        rung = list(pending)
        results: dict[int, IndexingResult] = {}
        positions: dict[int, int] = {}
        try:
            _run_program(command, working_directory, "dials.index")
        except ValueError as e:
            xia2_logger.debug(f"Indexing with method={method} failed:\n{e}")
        else:
            results, positions = _assess_rung(
                working_directory, tag, method, params, rung, n_strong
            )
        _record_log(f"index ({tag})", working_directory / f"dials.index.{tag}.log")

        for exposure in rung:
            result = results.get(exposure)
            if result is None:
                result = IndexingResult(
                    method=method,
                    n_strong=n_strong[exposure],
                    reason="no solution found",
                )
            attempts[exposure].append(result)
            xia2_logger.info(f"Exposure {exposure + 1}: {result.summary()}")
            if not result.accepted:
                continue
            accepted[exposure] = result
            locations[(exposure, method)] = (tag, positions[exposure])
            good_enough = (
                result.rmsd_px is None
                or params.target_rmsd_px is None
                or result.rmsd_px <= params.target_rmsd_px
            )
            if good_enough:
                # Good enough that a later rung would only cost time.
                pending.remove(exposure)
            # An exposure whose solution was merely acceptable stays pending, so
            # that a later rung can try to better it. Its accepted result is
            # kept either way, and the best of them wins below.

    if not accepted:
        return attempts

    # A later rung is walked because an earlier one was merely acceptable, not
    # because it is bound to be better, so each exposure keeps its best attempt
    # by RMSD rather than its last.
    chosen: dict[int, tuple[str, int]] = {}
    for exposure in sorted(accepted):
        solutions = [a for a in attempts[exposure] if a.accepted]
        best = min(solutions, key=lambda a: (a.rmsd_px is None, a.rmsd_px))
        if len(solutions) > 1 and attempts[exposure][-1] is not best:
            xia2_logger.info(
                f"Exposure {exposure + 1}: keeping the {best.method} solution, by RMSD"
            )
        # The accepted solution is reported last, as the outcome.
        attempts[exposure].remove(best)
        attempts[exposure].append(best)
        chosen[exposure] = locations[(exposure, best.method)]

    _merge_indexed(working_directory, chosen)
    return attempts


# How dials.refine parameterises the detector, always. Hierarchical rather
# than automatic so that the panels are free to move relative to one another,
# which is what a multi-panel TOF instrument needs: with one rigid detector,
# refinement absorbs a per-panel error into the crystal orientation instead of
# correcting it, and leaves residuals no mask or profile model can recover
# from. The level matters as much as the choice - on a flat hierarchy, which is
# what the three NMX panels hang off, level 0 *is* the whole detector as one
# body and refines exactly what automatic does (16 parameters, identical
# residuals); level 1 reaches the panels themselves (28 parameters).
DETECTOR_PANELS = "hierarchical"
DETECTOR_HIERARCHY_LEVEL = 1


def refine(working_directory: pathlib.Path) -> None:
    """Run dials.refine on the indexing solution, for every exposure at once."""
    command = [
        _executable("dials.refine"),
        "indexed.expt",
        "indexed.refl",
        f"refinement.parameterisation.detector.panels={DETECTOR_PANELS}",
        "refinement.parameterisation.detector.hierarchy_level"
        f"={DETECTOR_HIERARCHY_LEVEL}",
    ]
    xia2_logger.notice(banner("Refining"))  # type: ignore
    _run_program(command, working_directory, "dials.refine")
    _record_log("refine", working_directory / "dials.refine.log")


def _tof_integrate_command(
    params: IntegrationParams,
    method: str,
    mask: str,
    tag: str,
    experiments: str,
    reflections: str,
    pass_: IntegrationPass,
) -> list[str]:
    """The dials.tof_integrate command line for one method/mask combination."""
    command = [
        _executable("dials.tof_integrate"),
        experiments,
        reflections,
        f"method={method}",
        f"mask={mask}",
        f"integration_type={pass_.integration_type}",
        f"background_model={params.background_model}",
        f"ellipse_mask.scale={params.ellipse_mask_scale}",
        f"bbox_tof_padding={params.bbox_tof_padding}",
        f"bbox_xy_padding={params.bbox_xy_padding}",
        # Always: nothing after integration can apply the Lorentz factor,
        # and lawless has no LORENTZ keyword
        "corrections.lorentz=True",
        f"mp.nproc={params.nproc}",
        f"output.experiments=integrated_{tag}.expt",
        f"output.reflections=integrated_{tag}.refl",
        f"output.log=tof_integrate_{tag}.log",
    ]
    if params.phil:
        command.insert(1, os.fspath(params.phil))
    if pass_.d_min:
        command.append(f"calculated.dmin={pass_.d_min}")
    if pass_.wavelength_range:
        low, high = pass_.wavelength_range
        command.append(f"wavelength_range={low},{high}")
    if params.incident_run:
        command.append(f"corrections.incident_run={os.fspath(params.incident_run)}")
    if params.empty_run:
        command.append(f"corrections.empty_run={os.fspath(params.empty_run)}")
    if params.absorption:
        for name, value in params.absorption.items():
            command.append(f"corrections.absorption.target_spectrum.{name}={value}")
    return command


def _integration_quality(reflections: pathlib.Path) -> tuple[int, float]:
    """
    How well an integration run went: how many reflections were integrated by
    summation, and their mean I/sigma.

    Summation is used for the comparison whatever the method, because it is the
    one set of intensities that every run produces.
    """
    table = flex.reflection_table.from_file(reflections)
    table = table.select(table.get_flags(table.flags.integrated_sum))
    if not table.size():
        return 0, 0.0
    variance = table["intensity.sum.variance"]
    table = table.select(variance > 0)
    if not table.size():
        return 0, 0.0
    i_over_sigma = table["intensity.sum.value"] / flex.sqrt(
        table["intensity.sum.variance"]
    )
    return table.size(), float(flex.mean(i_over_sigma))


def _profile_fraction(reflections: pathlib.Path) -> float:
    """
    The fraction of integrated reflections that a profile fit succeeded for.

    A profile run can finish cleanly having fitted almost nothing, leaving an
    intensity.prf column that covers a handful of reflections, so success of
    the program is not success of the method.
    """
    table = flex.reflection_table.from_file(reflections)
    if "intensity.prf.value" not in table:
        return 0.0
    n_summed = table.get_flags(table.flags.integrated_sum).count(True)
    if not n_summed:
        return 0.0
    return table.get_flags(table.flags.integrated_prf).count(True) / n_summed


def _run_tof_integrate(
    working_directory: pathlib.Path,
    params: IntegrationParams,
    method: str,
    mask: str,
    experiments: str,
    reflections: str,
    pass_: IntegrationPass,
) -> tuple[str, pathlib.Path]:
    """
    Integrate with one mask, falling back to summation if a profile method
    fails outright.

    Profile fitting is not a per-reflection risk only: a single shoebox that
    the fit cannot bound aborts dials.tof_integrate, so the whole run is lost.
    Returns the method that actually produced the output.
    """
    tag = f"{method}_{mask}"
    if pass_.integration_type == "calculated":
        tag = f"{tag}_calculated"
    try:
        _run_program(
            _tof_integrate_command(
                params, method, mask, tag, experiments, reflections, pass_
            ),
            working_directory,
            "dials.tof_integrate",
        )
    except ValueError as e:
        fallback = params.fallback_method
        if method == fallback or fallback == NO_FALLBACK:
            raise
        xia2_logger.warning(
            f"Integration with method={method} failed, falling back to"
            f" {fallback}. The failure was:\n{str(e).strip().splitlines()[-1]}"
        )
        _record_log(
            f"{working_directory.name} integrate ({tag}, failed)",
            working_directory / f"tof_integrate_{tag}.log",
        )
        return _run_tof_integrate(
            working_directory, params, fallback, mask, experiments, reflections, pass_
        )

    _record_log(
        f"{working_directory.name} integrate ({tag})",
        working_directory / f"tof_integrate_{tag}.log",
    )
    integrated = working_directory / f"integrated_{tag}.refl"
    if not integrated.is_file():
        raise ValueError(f"dials.tof_integrate did not write {integrated.name}")

    if method != "summation":
        fraction = _profile_fraction(integrated)
        if fraction < params.min_profile_fraction:
            # Nothing to rerun: the summation intensities in this file are good
            # and the profile ones are kept as they are, so that both reach the
            # unmerged MTZ. What changes is which of them is exported as the
            # single SHELX intensity.
            xia2_logger.warning(
                f"{method} fitted only {fraction:.1%} of the reflections, fewer"
                f" than integration.min_profile_fraction"
                f" ({params.min_profile_fraction:.0%}), so the summation"
                " intensities will be the ones used"
            )
        else:
            xia2_logger.info(f"{method} fitted {fraction:.1%} of the reflections")
    return method, integrated


def _record_lorentz(reflections: pathlib.Path, applied: bool) -> None:
    """
    Record whether the Lorentz correction was applied, per reflection.

    dials.tof_integrate applies the correction to the intensities themselves
    and leaves no trace of it in the output beyond a line in its log, so the
    file cannot be asked later whether it has been corrected. The Lorentz
    factor must be applied exactly once, and the scaling programs downstream
    have their own keyword for it, so carry the answer with the data.
    """
    table = flex.reflection_table.from_file(reflections)
    table["lorentz_applied"] = flex.bool(table.size(), applied)
    table.as_file(reflections)


def _predicted_count(experiments, d_min: float) -> int:
    """
    How many reflections the geometry predicts out to d_min.

    Prediction is cheap next to integration - seconds against the hours that a
    hundred thousand shoeboxes would take - so the count is worth having before
    committing to a pass. It is an upper bound on what will be integrated,
    since dials.tof_integrate then drops overlapping reflections and any that
    fall outside the wavelength or time-of-flight range.
    """
    from dials.algorithms.spot_prediction import TOFReflectionPredictor

    total = 0
    for experiment in experiments:
        predictor = TOFReflectionPredictor(experiment, d_min)
        total += len(predictor.for_ub(experiment.crystal.get_A()))
    return total


def _d_min_for_count(experiments, target: int, d_min: float, predicted: int) -> float:
    """
    The resolution at which the geometry predicts about target reflections.

    Found by prediction rather than from the d^-3 that the number of
    reflections in a sphere follows, because a Laue-TOF experiment does not
    record a sphere: the wavelength band and the detector coverage both cut
    into it, so the analytic estimate comes out optimistic and a user who
    followed it would be refused a second time. It is only ever computed to
    put a number in that refusal, and prediction is a second or two.
    """
    low, high = d_min, d_min * (predicted / target) ** (1 / 3)
    for _ in range(8):
        if _predicted_count(experiments, high) <= target:
            break
        low, high = high, high * 1.2
    else:
        return high
    for _ in range(5):
        middle = 0.5 * (low + high)
        if _predicted_count(experiments, middle) <= target:
            high = middle
        else:
            low = middle
    return high


def _integrate_calculated(
    working_directory: pathlib.Path,
    params: IntegrationParams,
    method: str,
    mask: str,
    experiments: str,
    reflections: str,
    observed: pathlib.Path,
) -> tuple[str, pathlib.Path]:
    """
    Integrate every reflection the geometry predicts, not only those observed.

    This is the second pass of a calculated run. The first has already fixed
    the method and the mask on the spots that were found, which is cheap; this
    one measures everything that can have been recorded, most of it weak. That
    is the point: where to cut is a merging decision, taken on CC-half after
    scaling, and re-merging costs nothing where re-integrating costs
    everything.
    """
    expts = load.experiment_list(working_directory / experiments, check_format=False)
    limit, _ = geometric_d_min(expts[0])
    d_min = params.calculated_d_min or limit
    if params.calculated_d_min:
        xia2_logger.info(
            f"Integrating every predicted reflection to"
            f" integration.calculated_d_min = {d_min:.3f} A"
            + (
                f", beyond the {limit:.3f} A the geometry allows"
                if d_min < limit
                else f" ({limit:.3f} A is the geometric limit)"
            )
        )
    else:
        xia2_logger.info(
            f"Integrating every predicted reflection to the geometric limit of"
            f" {d_min:.3f} A"
        )

    # The band the data actually cover. Predictions outside it cannot have been
    # recorded, and integrating them adds noise with a Miller index on it. Only
    # the calculated pass needs this: an observed reflection is inside the band
    # by construction, having been found in the data.
    band = params.wavelength_range or wavelength_range(expts[0])

    n_predicted = _predicted_count(expts, d_min)
    n_observed = flex.reflection_table.from_file(observed).size()
    xia2_logger.info(
        f"{n_predicted} reflections predicted to {d_min:.3f} A, against"
        f" {n_observed} integrated from the spot list"
        f" ({n_predicted / max(n_observed, 1):.0f}x), over wavelengths"
        f" {band[0]:.2f}-{band[1]:.2f} A"
    )
    if n_predicted > params.max_predicted:
        suggestion = _d_min_for_count(expts, params.max_predicted, d_min, n_predicted)
        raise ValueError(
            f"{n_predicted} predicted reflections is more than"
            f" integration.max_predicted ({params.max_predicted}). Integrating"
            f" them would take hours. Either raise the guard, or set"
            f" integration.calculated_d_min to about {suggestion:.2f} A, which"
            f" is where this crystal reaches that many reflections."
        )

    method, integrated = _run_tof_integrate(
        working_directory,
        params,
        method,
        mask,
        experiments,
        reflections,
        IntegrationPass("calculated", d_min, band),
    )
    n_integrated = flex.reflection_table.from_file(integrated).size()
    xia2_logger.info(
        f"Integrated {n_integrated} of the {n_predicted} predicted reflections;"
        f" the {n_observed} from the spot list are kept in {observed.name}"
    )
    return method, integrated


def integrate(
    working_directory: pathlib.Path,
    params: IntegrationParams,
    choices: IntegrationChoices | None = None,
    experiments: str = "refined.expt",
    reflections: str = "refined.refl",
) -> pathlib.Path:
    """
    Run dials.tof_integrate, returning the integrated reflection file.

    With mask=both the data are integrated with each foreground mask and the
    better of the two is kept; that choice, and any fallback from a failed
    profile method, are recorded in choices and reused for the remaining
    orientations, so that every batch is integrated the same way.
    """
    choices = choices if choices is not None else IntegrationChoices()
    method = choices.method or params.method
    if choices.mask:
        masks = [choices.mask]
    elif params.mask == "both":
        masks = list(MASKS)
    else:
        masks = [params.mask]

    xia2_logger.notice(banner("Integrating"))  # type: ignore

    integrated_files = {}
    for mask in masks:
        method, integrated = _run_tof_integrate(
            working_directory,
            params,
            method,
            mask,
            experiments,
            reflections,
            IntegrationPass(),
        )
        integrated_files[mask] = integrated
    choices.method = method

    if len(integrated_files) > 1:
        quality = {
            mask: _integration_quality(f) for mask, f in integrated_files.items()
        }
        for mask, (n_integrated, i_over_sigma) in quality.items():
            xia2_logger.info(
                f"{mask}: {n_integrated} reflections integrated,"
                f" <I/sigma> {i_over_sigma:.1f}"
            )
        index = 1 if params.mask_metric == "i_over_sigma" else 0
        mask = max(quality, key=lambda m: quality[m][index])
        xia2_logger.info(f"Keeping the {mask} mask, by {params.mask_metric}")
    else:
        mask = masks[0]
    choices.mask = mask

    chosen = integrated_files[mask]
    if params.integration_type == "calculated":
        # The first pass is always over the observed spots: it is what chooses
        # the method and the mask, at a hundredth of the cost of doing that
        # comparison on every predicted reflection. Its output is left on disk
        # as a checkpoint, being a strict subset of what follows.
        # A profile method that survived the observed pass can still fail here,
        # where most of the shoeboxes hold no peak, so the method that came out
        # of this pass is the one the remaining orientations use.
        method, chosen = _integrate_calculated(
            working_directory, params, method, mask, experiments, reflections, chosen
        )
        choices.method = method

    integrated = working_directory / "integrated.refl"
    shutil.copyfile(chosen, integrated)
    shutil.copyfile(
        chosen.with_suffix(".expt"),
        working_directory / "integrated.expt",
    )
    _record_lorentz(integrated, True)
    FileHandler.record_data_file(str(working_directory / "integrated.expt"))
    FileHandler.record_data_file(str(integrated))
    return integrated


def combine(
    working_directory: pathlib.Path,
    integrated_directory: pathlib.Path,
    results: list[ExposureResult],
) -> tuple[pathlib.Path, pathlib.Path]:
    """
    Split the run's integrated output per exposure and write it ready to export.

    The exposures were integrated together, so this no longer concatenates
    files; it gives each exposure its own experiment identifier and its own
    batch, so that the exported intensities can be scaled per setting
    downstream.
    """
    indexed = [result for result in results if result.indexed]
    if not indexed:
        raise ValueError("No exposure was integrated, so there is nothing to combine")
    working_directory.mkdir(parents=True, exist_ok=True)
    xia2_logger.notice(banner("Combining"))  # type: ignore

    expts = load.experiment_list(
        integrated_directory / "integrated.expt", check_format=False
    )
    table = flex.reflection_table.from_file(integrated_directory / "integrated.refl")
    tables = parse_multiple_datasets([table])
    experiments = ExperimentList()
    experiments.extend(expts)
    # Which exposure each of the experiments came from. indexed.expt was built
    # in input order, and nothing since has reordered it.
    origins = [result.name for result in indexed][: len(tables)]
    experiments, tables = assign_unique_identifiers(experiments, tables)
    _match_intensity_columns(tables, origins)

    # dials.export numbers the Laue batches by imageset_id. Renumber it to the
    # position of the experiment in the list, so that each exposure exports as
    # its own batch rather than every one of them as batch 0.
    for i, (table, origin) in enumerate(zip(tables, origins)):
        table["imageset_id"] = flex.int(table.size(), i)
        xia2_logger.info(f"Batch {i}: {origin}, {table.size()} reflections")

    if len(tables) > 1:
        combined = flex.reflection_table.concat(tables)
    else:
        combined = tables[0]
    _report_lorentz(combined)

    experiments_file = working_directory / "combined.expt"
    reflections_file = working_directory / "combined.refl"
    experiments.as_file(experiments_file)
    combined.as_file(reflections_file)
    xia2_logger.info(
        f"Combined {len(experiments)} experiment(s), {combined.size()} reflections,"
        f" into {experiments_file.name} and {reflections_file.name}"
    )
    FileHandler.record_data_file(str(experiments_file))
    FileHandler.record_data_file(str(reflections_file))
    return experiments_file, reflections_file


def _intensity_choice(reflections: pathlib.Path, params: ExportParams) -> str:
    """
    Which intensities to export.

    dials.export refuses its own auto choice for SHELX when both summation and
    profile-fitted intensities are present, which is exactly what a successful
    profile run produces, so the choice is always made here.
    """
    if params.intensity != "auto":
        return params.intensity
    fraction = _profile_fraction(reflections)
    if fraction >= params.min_profile_fraction:
        return "profile"
    if fraction:
        xia2_logger.info(
            f"Profile-fitted intensities cover only {fraction:.1%} of the"
            " reflections, so the summation intensities are exported"
        )
    return "sum"


def _match_intensity_columns(
    tables: list[flex.reflection_table], origins: list[str]
) -> None:
    """
    Make every table carry the same intensity columns, so they can be combined.

    A profile method that aborted on one orientation and worked on another
    leaves tables that cannot be concatenated, and an MTZ cannot hold a column
    for only some of its observations. The profile-fitted intensities are the
    ones that go, since summation is always there.
    """
    with_profile = [
        origin
        for table, origin in zip(tables, origins)
        if "intensity.prf.value" in table
    ]
    if not with_profile or len(with_profile) == len(tables):
        return
    missing = [
        origin
        for table, origin in zip(tables, origins)
        if "intensity.prf.value" not in table
    ]
    xia2_logger.warning(
        f"{', '.join(missing)} {'has' if len(missing) == 1 else 'have'} no"
        " profile-fitted intensities, so they are dropped from"
        f" {', '.join(with_profile)} too and only summation intensities are"
        " combined."
    )
    for table in tables:
        for column in ("intensity.prf.value", "intensity.prf.variance"):
            if column in table:
                del table[column]


def _lorentz_applied(reflections_file: pathlib.Path) -> bool:
    """
    Whether the data already carries the Lorentz correction.

    This pipeline always applies it at integration, so the answer is yes for
    anything it produced; the check exists for files it did not, which it can
    be pointed at by rerunning in an old directory. Without the column the
    answer is no, which is what dials.tof_integrate defaults to and so what
    data integrated elsewhere most likely carries. The combine step has already
    warned about the missing provenance.
    """
    table = flex.reflection_table.from_file(reflections_file)
    if "lorentz_applied" not in table:
        return False
    return bool(set(table["lorentz_applied"]).pop())


def _report_lorentz(table: flex.reflection_table) -> None:
    """
    Say whether the combined data carry the Lorentz correction.

    The correction has to be applied exactly once, here or by the scaling
    program, so the answer decides what is passed downstream. Orientations that
    disagree cannot be scaled together.
    """
    if "lorentz_applied" not in table:
        xia2_logger.warning(
            "The combined data do not record whether the Lorentz correction was"
            " applied. It must be applied exactly once - check before scaling."
        )
        return
    applied = set(table["lorentz_applied"])
    if len(applied) > 1:
        raise ValueError(
            "Some orientations were integrated with the Lorentz correction and"
            " some without, so they cannot be scaled together. Reintegrate them"
            " the same way."
        )
    if applied.pop():
        xia2_logger.info(
            "The Lorentz correction was applied during integration, which is"
            " the only place this pipeline can take it."
        )
    else:
        xia2_logger.warning(
            "The Lorentz correction was not applied during integration, and no"
            " step after this one can apply it. These data were not integrated"
            " by this pipeline, which always applies it; reintegrate them."
        )


def export_shelx(
    working_directory: pathlib.Path,
    params: ExportParams,
    experiments: str = "combined.expt",
    reflections: str = "combined.refl",
) -> pathlib.Path:
    """Run dials.export format=shelx, returning the exported hkl file."""
    intensity = _intensity_choice(working_directory / reflections, params)
    if not _lorentz_applied(working_directory / reflections):
        xia2_logger.warning(
            "The exported intensities do not carry the Lorentz correction, and"
            " nothing else applies it - lawless has no LORENTZ keyword. These"
            " data were not integrated by this pipeline; reintegrate them."
        )
    command = [
        _executable("dials.export"),
        experiments,
        reflections,
        "format=shelx",
        f"intensity={intensity}",
        f"shelx.composition={params.composition}",
    ]
    if params.phil:
        command.insert(1, os.fspath(params.phil))

    xia2_logger.notice(banner("Exporting"))  # type: ignore
    _run_program(command, working_directory, "dials.export")
    _record_log("export", working_directory / "dials.export.log")

    hklout = working_directory / "dials.hkl"
    if not hklout.is_file():
        raise ValueError("dials.export did not write dials.hkl")
    xia2_logger.info(f"Exported the {intensity} intensities to {hklout.name}")
    FileHandler.record_data_file(str(hklout))
    ins = working_directory / "dials.ins"
    if ins.is_file():
        FileHandler.record_data_file(str(ins))
    return hklout


def unmerged_mtz(
    working_directory: pathlib.Path,
    experiments_file: pathlib.Path,
    reflections_file: pathlib.Path,
) -> pathlib.Path:
    """
    Write the combined data as an unmerged Laue MTZ, for pointless and lawless.

    dials.export cannot do this - format=mtz is refused for time-of-flight data
    - so xia2 writes the file itself.
    """
    from xia2.Modules.Laue_TOF.laue_tof_mtz import write_unmerged_mtz

    xia2_logger.notice(banner("Writing the unmerged MTZ"))  # type: ignore
    mtz_file = write_unmerged_mtz(
        experiments_file, reflections_file, working_directory / "unmerged.mtz"
    )
    FileHandler.record_data_file(str(mtz_file))
    return mtz_file


def _report(results: list[ExposureResult]) -> None:
    """Summarise what happened to each exposure."""
    xia2_logger.notice(banner("Summary"))  # type: ignore
    for result in results:
        if result.integrated:
            state = "integrated"
        elif result.indexed:
            state = "indexed, not integrated"
        elif result.n_strong:
            state = "not indexed"
        else:
            state = "no strong spots"
        line = f"{result.name}: {state}"
        if result.indexing and result.indexing.accepted:
            line += f" ({result.indexing.summary()})"
        xia2_logger.info(line)


def run_data_integration(
    root_working_directory: pathlib.Path, setup: LaueTOFSetup
) -> list[ExposureResult]:
    """
    Process one run: every exposure of it, through the DIALS steps together.

    The exposures are stationary measurements of one crystal at different
    settings, so they belong to one run and are handed to each DIALS program as
    a single multi-experiment set. Indexing is the step that has to know they
    are separate measurements rather than one sweep, which is what
    joint_indexing=False tells it.
    """
    inputs = [("image", image) for image in setup.file_input.images]
    inputs += [("template", t) for t in setup.file_input.templates]
    inputs += [("directory", d) for d in setup.file_input.directories]

    run_directory = root_working_directory / "run"
    run_directory.mkdir(parents=True, exist_ok=True)
    xia2_logger.notice(banner(f"Processing {len(inputs)} exposure(s)"))  # type: ignore
    for i, (_, path) in enumerate(inputs, start=1):
        xia2_logger.info(f"Exposure {i}: {path}")

    results = [
        ExposureResult(name=f"exposure_{i}", image=path)
        for i, (_, path) in enumerate(inputs, start=1)
    ]

    run_import(run_directory, inputs, setup.file_input)

    if "find_spots" not in setup.options.steps:
        return results
    counts = find_spots(run_directory, setup.spotfinding_params, len(inputs))
    for result, count in zip(results, counts):
        result.n_strong = count

    if "index" not in setup.options.steps:
        return results
    attempts = index(run_directory, setup.indexing_params, counts)
    for result, exposure_attempts in zip(results, attempts):
        result.attempts = exposure_attempts
        result.indexing = exposure_attempts[-1] if exposure_attempts else None
        if not result.indexed:
            reason = result.indexing.reason if result.indexing else "no attempt made"
            xia2_logger.warning(f"{result.name} did not index: {reason}")

    indexed = [result for result in results if result.indexed]
    if not indexed:
        xia2_logger.warning("No exposure indexed, so there is nothing to refine")
        _report(results)
        return results
    if len(indexed) < len(results):
        xia2_logger.warning(
            f"{len(results) - len(indexed)} of {len(results)} exposures did not"
            " index and are left out of the rest of the run"
        )

    if "refine" not in setup.options.steps:
        _report(results)
        return results
    refine(run_directory)

    if "integrate" not in setup.options.steps:
        _report(results)
        return results
    integrate(run_directory, setup.integration_params)
    for result in indexed:
        result.integrated = True

    _report(results)

    if "combine" not in setup.options.steps:
        return results

    scale_directory = root_working_directory / "scale"
    experiments_file, reflections_file = combine(
        scale_directory, run_directory, results
    )
    if "export" in setup.options.steps:
        export_shelx(scale_directory, setup.export_params)
    if "unmerged_mtz" not in setup.options.steps:
        return results
    mtz_file = unmerged_mtz(scale_directory, experiments_file, reflections_file)

    if "pointless" in setup.options.steps or "lawless" in setup.options.steps:
        from xia2.Modules.Laue_TOF.laue_tof_scale import scale

        scale(
            scale_directory,
            mtz_file,
            setup.pointless_params,
            setup.lawless_params,
            _lorentz_applied(reflections_file),
            setup.options.steps,
            report_directory=root_working_directory,
        )
    return results
