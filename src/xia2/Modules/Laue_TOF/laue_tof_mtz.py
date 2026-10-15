"""
Check the unmerged Laue MTZ that dials.export writes.

The plan's "The unmerged MTZ" table lists what fails silently downstream;
validate_unmerged_mtz checks each of them.
"""

from __future__ import annotations

import pathlib

import gemmi
import numpy

# MTZ batch header word offsets, as used by CCP4 and gemmi
I_LDTYPE = 14  # 1 = 2D oscillation, 2 = 3D area detector, 3 = Laue
F_UMAT = 6  # 6..14
F_BSCALE = 43


def check_wavelength_column(mtz_file: pathlib.Path) -> str:
    """
    The name of the per-observation wavelength column in an MTZ, or "".

    lawless stops on LAUE without it, and takes the first of these labels it finds.
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
                "M/ISYM is 0 for some observations, so the file reads as merged"
            )

    if not mtz.batches:
        problems.append("no batch records, so lawless cannot tell the exposures apart")
    else:
        if expected_batches is not None and len(mtz.batches) != expected_batches:
            problems.append(
                f"{len(mtz.batches)} batch records for {expected_batches} orientations"
            )
        for batch in mtz.batches:
            if batch.ints[I_LDTYPE] != 3:
                problems.append(f"batch {batch.number} is not LDTYPE 3 (Laue)")
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
