"""
Finish and check the unmerged Laue MTZ that dials.export writes.

dials.export writes one batch per exposure and a LAMBDA column, but leaves four
things pointless and lawless need: the Laue batch type, a non-zero phi range, each
exposure's own wavelength band, and intensities on a counts-like scale. The plan's "The
unmerged MTZ" table lists what fails silently, and validate_unmerged_mtz checks
each of them.
"""

from __future__ import annotations

import logging
import math
import pathlib

import gemmi
import numpy

xia2_logger = logging.getLogger(__name__)

# MTZ batch header word offsets, as used by CCP4 and gemmi
I_LDTYPE = 14  # 1 = 2D oscillation, 2 = 3D area detector, 3 = Laue
F_UMAT = 6  # 6..14
F_PHISTT = 36
F_PHIEND = 37
F_BSCALE = 43
F_PHIRANGE = 47
F_ALAMBD = 86
F_DELAMB = 87

# What the intensities are scaled towards; see _counts_like_scale
COUNTS_LIKE_MEDIAN = 1.0e3

INTENSITY_COLUMNS = ("I", "SIGI", "IPR", "SIGIPR")


def _counts_like_scale(intensity: numpy.ndarray) -> float:
    """
    A power of ten that puts the intensities on a counts-like scale.

    The Lorentz correction leaves them of order 1e-3, where pointless fails its
    normalisation outright ("NormaliseInRuns::apply Unused run 1").
    """
    finite = intensity[numpy.isfinite(intensity)]
    finite = numpy.abs(finite[finite != 0])
    if not finite.size:
        return 1.0
    return float(10.0 ** round(math.log10(COUNTS_LIKE_MEDIAN / numpy.median(finite))))


def finish_unmerged_mtz(mtz_file: pathlib.Path, expected_batches: int) -> None:
    """Mark the batches as Laue, give each a phi range and its band, and rescale."""
    mtz = gemmi.read_mtz_file(str(mtz_file))
    data = numpy.array(mtz, copy=True)
    labels = [column.label for column in mtz.columns]

    scale = _counts_like_scale(data[:, labels.index("I")])
    if scale != 1.0:
        xia2_logger.info(
            f"Scaling the intensities by {scale:g} to put them on the counts-like"
            " scale pointless and lawless expect; a constant factor changes no"
            " statistic."
        )
        for label in INTENSITY_COLUMNS:
            if label in labels:
                data[:, labels.index(label)] *= scale
        mtz.set_data(data)

    batch_of_row = data[:, labels.index("BATCH")]
    wavelength = data[:, labels.index("LAMBDA")]
    for batch in mtz.batches:
        batch.ints[I_LDTYPE] = 3
        # lawless drops SECONDARY and ABSORPTION, with only a warning, on PHIRANGE 0
        batch.floats[F_PHISTT] = float(batch.number)
        batch.floats[F_PHIEND] = float(batch.number) + 1.0
        batch.floats[F_PHIRANGE] = 1.0
        band = wavelength[batch_of_row == batch.number]
        if band.size:
            batch.floats[F_ALAMBD] = float(band.mean())
            batch.floats[F_DELAMB] = 0.5 * float(band.max() - band.min())
        else:
            xia2_logger.warning(f"Batch {batch.number} has no reflections")

    mtz.write_to_file(str(mtz_file))
    validate_unmerged_mtz(mtz_file, expected_batches)


def check_wavelength_column(mtz_file: pathlib.Path) -> str:
    """
    The name of the per-observation wavelength column in an MTZ, or "".

    pointless builds up to 1.13.6 drop it, which makes LAUE fatal in lawless.
    lawless takes the first of these labels it finds.
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    labels = {column.label.upper(): column.label for column in mtz.columns}
    for candidate in ("LAMBDA", "LAM", "WAVELENGTH"):
        if candidate in labels:
            return labels[candidate]
    return ""


def validate_unmerged_mtz(
    mtz_file: pathlib.Path, expected_batches: int | None = None
) -> None:
    """
    Check the things that fail silently downstream.

    Raises ValueError naming the specific cause, because the symptoms - a
    collapsed merged file, one batch reported for everything, a flat secondary
    correction - do not point back here.
    """
    mtz = gemmi.read_mtz_file(str(mtz_file))
    problems = []

    if not check_wavelength_column(mtz_file):
        problems.append("no LAMBDA column, so lawless cannot run LAUE")

    labels = [column.label for column in mtz.columns]
    for required in ("I", "SIGI", "M/ISYM", "BATCH"):
        if required not in labels:
            problems.append(f"no {required} column")

    if "M/ISYM" in labels:
        misym = numpy.asarray(mtz.column_with_label("M/ISYM").array)
        if misym.size and not bool((misym > 0).all()):
            problems.append(
                "M/ISYM is 0 for some observations, so pointless will read the"
                " file as merged and collapse them"
            )

    if not mtz.batches:
        problems.append(
            "no batch records, so pointless will report one batch whatever the"
            " BATCH column says"
        )
    else:
        if expected_batches is not None and len(mtz.batches) != expected_batches:
            problems.append(
                f"{len(mtz.batches)} batch records for {expected_batches} orientations"
            )
        for batch in mtz.batches:
            if batch.ints[I_LDTYPE] != 3:
                problems.append(f"batch {batch.number} is not LDTYPE 3 (Laue)")
            if not batch.floats[F_PHIRANGE]:
                problems.append(
                    f"batch {batch.number} has PHIRANGE 0, so lawless drops"
                    " SECONDARY and ABSORPTION"
                )
            if batch.floats[F_BSCALE] != 1.0:
                problems.append(f"batch {batch.number} has BSCALE != 1")
            umat = numpy.array(
                [batch.floats[F_UMAT + i] for i in range(9)], dtype=float
            ).reshape(3, 3)
            if numpy.allclose(umat, numpy.eye(3)):
                problems.append(
                    f"batch {batch.number} has an identity UMAT, which makes"
                    " SCALES SECONDARY and ABSORPTION refine flat"
                )
            if numpy.linalg.det(umat) <= 0.9:
                problems.append(
                    f"batch {batch.number} has an invalid UMAT (det"
                    f" {numpy.linalg.det(umat):.3f})"
                )

    if problems:
        raise ValueError(
            f"{mtz_file} is not a usable unmerged Laue MTZ:\n  " + "\n  ".join(problems)
        )
