from __future__ import annotations

import errno
import functools
import importlib
import logging
import os
import pathlib
from dataclasses import dataclass, field

import iotbx.phil
from cctbx import crystal, sgtbx, uctbx
from dials.util.system import CPU_COUNT
from libtbx import Auto

from xia2.Modules.Laue_TOF.laue_tof_phil import *
from xia2.Modules.Laue_TOF.util import report_timing

xia2_logger = logging.getLogger(__name__)


# What a failed profile method can fall back to, besides another method
NO_FALLBACK = "none"


def _resolved_file(value: str | None) -> pathlib.Path | None:
    """Resolve a phil path option, checking that the file exists."""
    if not value:
        return None
    path = pathlib.Path(value).resolve()
    if not path.is_file():
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), os.fspath(path)
        )
    return path


@functools.cache
def _program_choices(module: str, path: str) -> tuple[str, ...]:
    """
    The values a choice parameter of another program's phil allows.

    Asked of the program rather than listed here, so that a DIALS which gains
    or loses a method needs no change on this side, and a name this pipeline
    would have accepted cannot turn out to be one the program has never heard
    of.
    """
    try:
        scope = importlib.import_module(module).phil_scope
    except (ImportError, AttributeError) as e:
        raise ValueError(f"Could not ask {module} which methods it offers: {e}")
    definition = scope.get(path).objects[0]
    # A choice marks its default with a leading asterisk
    return tuple(word.value.lstrip("*") for word in definition.words)


def indexing_methods() -> tuple[str, ...]:
    """The indexing methods dials.index offers."""
    return _program_choices("dials.command_line.index", "indexing.method")


def integration_methods() -> tuple[str, ...]:
    """The integration methods dials.tof_integrate offers."""
    return _program_choices("dials.command_line.tof_integrate", "method")


def _format_class(filename: pathlib.Path):
    """The dxtbx format class that claims a file, or None if none does."""
    from dxtbx.format.Registry import get_format_class_for_file

    try:
        return get_format_class_for_file(os.fspath(filename))
    except Exception:
        return None


def readable_images(filename: pathlib.Path) -> int:
    """
    How many images dxtbx can read from a file, or zero if it cannot read it.
    """
    format_class = _format_class(filename)
    if format_class is None:
        return 0
    try:
        return int(format_class(os.fspath(filename)).get_num_images())
    except Exception:
        return 0


def check_inputs_readable(file_input: FileInput) -> None:
    """
    Stop before any work if dials.import can not read an input file.
    """
    unreadable = []
    for image in file_input.images:
        path = pathlib.Path(image)
        if not path.is_file():
            raise FileNotFoundError(
                errno.ENOENT, os.strerror(errno.ENOENT), os.fspath(path)
            )
        if not readable_images(path):
            unreadable.append(image)
    if unreadable:
        raise ValueError(
            "dials.import cannot read "
            + ", ".join(unreadable)
            + ". Raw event data has to be reduced to time-of-flight bins"
            " before this pipeline can process it."
        )


@dataclass
class FileInput:
    images: list[str] = field(default_factory=list)
    templates: list[str] = field(default_factory=list)
    directories: list[str] = field(default_factory=list)
    vanadium_run: pathlib.Path | None = None
    empty_run: pathlib.Path | None = None
    mask: pathlib.Path | None = None
    detector_calibration: pathlib.Path | None = None
    import_phil: pathlib.Path | None = None

    def resolve_paths(self) -> None:
        for filetype in (self.images, self.templates):
            for i, obj in enumerate(filetype):
                filetype[i] = str(pathlib.Path(obj).resolve())

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> FileInput:
        file_input = cls()
        if params.image:
            file_input.images = list(params.image)
        elif params.template:
            file_input.templates = list(params.template)
        elif params.directory:
            file_input.directories = [
                str(pathlib.Path(i).resolve()) for i in params.directory
            ]
        else:
            raise ValueError(
                "No input data identified (use image=, template= or directory=)"
            )
        file_input.resolve_paths()
        file_input.vanadium_run = _resolved_file(params.vanadium_run)
        file_input.empty_run = _resolved_file(params.empty_run)
        file_input.mask = _resolved_file(params.mask)
        file_input.detector_calibration = _resolved_file(params.detector_calibration)
        file_input.import_phil = _resolved_file(params.dials_import.phil)
        if file_input.empty_run and not file_input.vanadium_run:
            raise ValueError(
                "An empty_run was given without a vanadium_run. The empty run is"
                " only used to correct the incident spectrum, so both are needed."
            )
        return file_input


@dataclass
class SpotfindingParams:
    min_spot_size: int = 15
    max_spot_size: int | None = None
    threshold_algorithm: str = "dispersion_extended"
    max_strong: int = 35000
    nproc: int = 1
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> SpotfindingParams:
        return cls(
            params.spotfinding.min_spot_size,
            params.spotfinding.max_spot_size,
            params.spotfinding.threshold_algorithm,
            params.spotfinding.max_strong,
            params.multiprocessing.nproc,
            _resolved_file(params.spotfinding.phil),
        )


@dataclass
class IndexingParams:
    space_group: sgtbx.space_group | None = None
    unit_cell: uctbx.unit_cell | None = None
    methods: list[str] = field(default_factory=lambda: ["fft3d"])
    max_cell: float | None = None
    d_min_start: float | None = 4.0
    min_spots: int = 10
    min_indexed: int = 50
    min_indexed_fraction: float = 0.25
    max_rmsd_px: float | None = 5.0
    target_rmsd_px: float | None = 2.0
    seed_from_first: bool = True
    outlier_algorithm: str = "mcd"
    outlier_iqr_multiplier: float = 0.8
    outlier_separate_panels: bool = True
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> IndexingParams:
        indexing = params.indexing
        if indexing.unit_cell and params.space_group:
            try:
                _ = crystal.symmetry(
                    unit_cell=indexing.unit_cell,
                    space_group_info=params.space_group,
                    assert_is_compatible_unit_cell=True,
                )
            except AssertionError as e:
                raise ValueError(e)

        if indexing.method:
            methods = [indexing.method]
        else:
            methods = list(indexing.ladder)
        if not methods:
            raise ValueError(
                "No indexing method given (set indexing.method or indexing.ladder)"
            )
        available = indexing_methods()
        for method in methods:
            if method not in available:
                raise ValueError(
                    f"dials.index has no indexing method {method}."
                    f" It offers {', '.join(available)}"
                )

        max_cell = indexing.max_cell
        if max_cell is None:
            # Always fix max_cell. The DIALS nearest-neighbour estimate is
            # unreliable on the sparse spot lists Laue-TOF data gives, and it
            # fails quietly: an estimate that comes out too small returns a
            # truncated cell which indexes fewer spots but refines well enough
            # to pass every acceptance test here.
            if indexing.unit_cell:
                max_cell = indexing.max_cell_multiplier * max(
                    indexing.unit_cell.parameters()[:3]
                )
            else:
                max_cell = indexing.max_cell_fallback

        return cls(
            params.space_group,
            indexing.unit_cell,
            methods,
            max_cell,
            indexing.d_min_start,
            indexing.min_spots,
            indexing.min_indexed,
            indexing.min_indexed_fraction,
            indexing.max_rmsd_px,
            indexing.target_rmsd_px,
            indexing.seed_from_first,
            indexing.outlier.algorithm,
            indexing.outlier.iqr_multiplier,
            indexing.outlier.separate_panels,
            _resolved_file(indexing.phil),
        )


@dataclass
class IntegrationParams:
    method: str = "profile_1d_ibix"
    fallback_method: str = "summation"
    min_profile_fraction: float = 0.5
    integration_type: str = "observed"
    calculated_d_min: float | None = None
    max_predicted: int = 250000
    background_model: str = "linear3d"
    mask: str = "both"
    mask_metric: str = "i_over_sigma"
    ellipse_mask_scale: float = 1.0
    wavelength_range: tuple[float, float] | None = None
    bbox_tof_padding: int = 2
    bbox_xy_padding: int = 1
    incident_run: pathlib.Path | None = None
    empty_run: pathlib.Path | None = None
    absorption: dict[str, float] | None = None
    nproc: int = 1
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(
        cls, params: iotbx.phil.scope_extract, file_input: FileInput
    ) -> IntegrationParams:
        integration = params.integration
        available = integration_methods()
        if integration.method not in available:
            raise ValueError(
                f"dials.tof_integrate has no integration method"
                f" {integration.method}. It offers {', '.join(available)}"
            )
        # phil reads the bare word none as None for a str parameter
        fallback_method = integration.fallback_method or NO_FALLBACK
        if fallback_method not in available + (NO_FALLBACK,):
            raise ValueError(
                f"integration.fallback_method {fallback_method} is"
                f" neither {NO_FALLBACK} nor one of the methods"
                f" dials.tof_integrate offers ({', '.join(available)})"
            )
        if (
            integration.calculated_d_min is not None
            and integration.calculated_d_min <= 0
        ):
            raise ValueError(
                "integration.calculated_d_min must be a positive resolution in"
                " Angstrom, or None to use the limit the geometry allows"
            )

        absorption = None
        if integration.absorption.enabled:
            absorption = {
                "sample_number_density": integration.absorption.sample_number_density,
                "sample_radius": integration.absorption.sample_radius,
                "scattering_x_section": integration.absorption.scattering_x_section,
                "absorption_x_section": integration.absorption.absorption_x_section,
            }
            missing = [name for name, value in absorption.items() if value is None]
            if missing:
                raise ValueError(
                    "integration.absorption.enabled=True requires the target"
                    f" spectrum parameters, but {', '.join(missing)} not set"
                )

        wavelength_range = None
        if integration.wavelength_range:
            wavelength_range = tuple(integration.wavelength_range)

        return cls(
            integration.method,
            fallback_method,
            integration.min_profile_fraction,
            integration.integration_type,
            integration.calculated_d_min,
            integration.max_predicted,
            integration.background_model,
            integration.mask,
            integration.mask_metric,
            integration.ellipse_mask_scale,
            wavelength_range,
            integration.bbox_tof_padding,
            integration.bbox_xy_padding,
            file_input.vanadium_run,
            file_input.empty_run,
            absorption,
            params.multiprocessing.nproc,
            _resolved_file(integration.phil),
        )


@dataclass
class ExportParams:
    intensity: str = "auto"
    min_profile_fraction: float = 0.5
    composition: str = "CH"
    phil: pathlib.Path | None = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> ExportParams:
        return cls(
            params.output.intensity,
            params.integration.min_profile_fraction,
            params.output.composition,
            _resolved_file(params.output.phil),
        )


@dataclass
class PointlessParams:
    executable: pathlib.Path | None = None
    lauegroup: str | None = None
    space_group: sgtbx.space_group | None = None
    keywords: list[str] = field(default_factory=list)

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> PointlessParams:
        pointless = params.scaling.pointless
        space_group = pointless.space_group or params.space_group
        lauegroup = pointless.lauegroup
        if lauegroup is None and space_group:
            group = (
                space_group.group() if hasattr(space_group, "group") else space_group
            )
            lauegroup = (
                sgtbx.space_group_info(
                    group=group.build_derived_reflection_intensity_group(
                        anomalous_flag=False
                    )
                )
                .type()
                .lookup_symbol()
            )
        return cls(
            _resolved_file(pointless.executable),
            lauegroup,
            space_group,
            [line for line in pointless.keywords if line],
        )


@dataclass
class LawlessParams:
    executable: pathlib.Path | None = None
    probe: str = "neutron"
    normalisation: str = "gpr"
    lam_min: float | None = None
    lam_max: float | None = None
    lam_ref: float | None = None
    chebyshev_degree: int = 6
    gpr_bins: int | None = None
    scales: str = "BATCH BFACTOR ON SECONDARY 0"
    sdcorrection: str = "NOREFINE 1.0 0.0 0.0"
    anomalous: bool = False
    resolution: tuple[float, float] | None = None
    lambda_only: bool = False
    keywords: list[str] = field(default_factory=list)

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> LawlessParams:
        lawless = params.scaling.lawless
        return cls(
            _resolved_file(lawless.executable),
            lawless.probe,
            lawless.normalisation,
            lawless.lam_min,
            lawless.lam_max,
            lawless.lam_ref,
            lawless.chebyshev_degree,
            lawless.gpr_bins,
            lawless.scales,
            lawless.sdcorrection,
            lawless.anomalous,
            tuple(lawless.resolution) if lawless.resolution else None,
            lawless.lambda_only,
            [line for line in lawless.keywords if line],
        )


@dataclass
class AlgorithmParams:
    steps: list[str] = field(default_factory=list)
    nproc: int = 1

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract) -> AlgorithmParams:
        return cls(list(params.workflow.steps), params.multiprocessing.nproc)


@dataclass
class LaueTOFSetup:
    """The fully resolved input of a run, as processed from the phil scope."""

    file_input: FileInput
    spotfinding_params: SpotfindingParams
    indexing_params: IndexingParams
    integration_params: IntegrationParams
    export_params: ExportParams
    pointless_params: PointlessParams
    lawless_params: LawlessParams
    options: AlgorithmParams


def _log_setup(setup: LaueTOFSetup) -> None:
    """Report what was resolved from the phil scope, before any work is done."""
    for image in setup.file_input.images:
        n_images = readable_images(pathlib.Path(image))
        xia2_logger.info(f"Input {image}: {n_images} image(s)")

    indexing = setup.indexing_params
    if len(indexing.methods) > 1:
        xia2_logger.info(
            "Indexing methods to try, in order: " + ", ".join(indexing.methods)
        )
    else:
        xia2_logger.info(f"Indexing method: {indexing.methods[0]}")
    if indexing.max_cell:
        if indexing.unit_cell:
            source = "from the longest edge of the given unit cell"
        else:
            source = (
                "a fixed default, since no unit cell was given. Raise"
                " indexing.max_cell, or give indexing.unit_cell, if the cell has"
                " a longer edge than this"
            )
        xia2_logger.info(f"Using max_cell = {indexing.max_cell:.1f} Angstrom, {source}")

    integration = setup.integration_params
    if integration.incident_run:
        xia2_logger.info(
            f"Normalising by the incident spectrum from {integration.incident_run}"
        )
        if not integration.empty_run:
            xia2_logger.warning(
                "No empty_run given, so the incident spectrum will not be"
                " background corrected."
            )
    else:
        xia2_logger.warning(
            "No vanadium_run given: the intensities will not be normalised by the"
            " incident spectrum."
        )
    xia2_logger.info(
        "The TOF Lorentz correction is applied during integration, per pixel,"
        " which is the only place this pipeline can take it"
    )
    if integration.mask == "both":
        xia2_logger.info(
            "Integrating with each foreground mask and keeping the better by"
            f" {integration.mask_metric}"
        )
    xia2_logger.info(
        f"Integration method: {integration.method}"
        + (
            ""
            if integration.method == "summation"
            or not integration.fallback_method
            or integration.fallback_method == NO_FALLBACK
            else f", falling back to {integration.fallback_method} if it fails"
        )
    )
    if integration.integration_type == "calculated":
        xia2_logger.info(
            "Every predicted reflection will be integrated, to "
            + (
                f"{integration.calculated_d_min} A"
                if integration.calculated_d_min
                else "the limit the geometry allows"
            )
            + ", after a first pass over the observed spots (at most"
            f" {integration.max_predicted} reflections, see"
            " integration.max_predicted)"
        )


def setup_from_phil(params: iotbx.phil.scope_extract) -> LaueTOFSetup:
    """Turn the phil scope into the set of parameter objects the steps take."""
    if params.multiprocessing.nproc is Auto:
        params.multiprocessing.nproc = CPU_COUNT

    file_input = FileInput.from_phil(params)
    check_inputs_readable(file_input)
    setup = LaueTOFSetup(
        file_input=file_input,
        spotfinding_params=SpotfindingParams.from_phil(params),
        indexing_params=IndexingParams.from_phil(params),
        integration_params=IntegrationParams.from_phil(params, file_input),
        export_params=ExportParams.from_phil(params),
        pointless_params=PointlessParams.from_phil(params),
        lawless_params=LawlessParams.from_phil(params),
        options=AlgorithmParams.from_phil(params),
    )
    _log_setup(setup)
    return setup


@report_timing
def run_xia2_laue_tof(
    root_working_directory: pathlib.Path, params: iotbx.phil.scope_extract
) -> LaueTOFSetup:
    """
    Run the Laue-TOF processing pipeline.

    Every exposure of the run goes through the DIALS steps together, in one
    working directory: import, find spots, index, refine and integrate. They
    are separate stationary measurements of one crystal, so indexing treats
    each on its own while the other steps see the whole set.
    """
    setup = setup_from_phil(params)

    try:
        from xia2.Modules.Laue_TOF.data_integration import run_data_integration
    except ImportError as e:
        raise NotImplementedError(
            "The processing steps are not implemented yet; this module currently"
            " only resolves and validates the input parameters."
        ) from e

    run_data_integration(root_working_directory, setup)
    return setup
