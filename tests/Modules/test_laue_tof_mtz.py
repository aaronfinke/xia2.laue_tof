from __future__ import annotations

import gemmi
import numpy
import pytest

from xia2.Modules.Laue_TOF.laue_tof_mtz import I_LDTYPE, validate_unmerged_mtz
from xia2.Modules.Laue_TOF.laue_tof_report import wavelength_band_plot


def _exported_mtz(path, ldtype=3, misym=1):
    """An MTZ shaped like dials.export's polychromatic output."""
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.SpaceGroup("P 1 21 1")
    mtz.cell = gemmi.UnitCell(43, 42, 73, 90, 104.5, 90)
    mtz.add_dataset("FROMDIALS")
    for label, kind in (("M/ISYM", "Y"), ("BATCH", "B")):
        mtz.add_column(label, kind, dataset_id=0)
    for label, kind in (("I", "J"), ("SIGI", "Q"), ("LAMBDA", "R")):
        mtz.add_column(label, kind, dataset_id=1)
    mtz.set_data(
        numpy.array(
            [
                [1, 0, 0, misym, 1, 2e-3, 1e-4, 2.0],
                [0, 1, 0, 2, 1, 4e-3, 2e-4, 3.0],
                [0, 0, 1, 1, 101, 3e-3, 1e-4, 3.5],
            ],
            dtype=numpy.float32,
        )
    )
    umat = numpy.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    for number in (1, 101):
        batch = gemmi.Mtz.Batch()
        batch.number = number
        batch.ints[I_LDTYPE] = ldtype
        batch.floats[43] = 1.0
        for i, value in enumerate(umat.T.reshape(9)):
            batch.floats[6 + i] = value
        batch.dataset_id = 1
        mtz.batches.append(batch)
    mtz.write_to_file(str(path))
    return path


def test_validate_unmerged_mtz(tmp_path):
    validate_unmerged_mtz(_exported_mtz(tmp_path / "ok.mtz"), expected_batches=2)
    with pytest.raises(ValueError, match="2 batch records for 3"):
        validate_unmerged_mtz(_exported_mtz(tmp_path / "ok.mtz"), expected_batches=3)
    with pytest.raises(ValueError, match="LDTYPE 3"):
        validate_unmerged_mtz(_exported_mtz(tmp_path / "ld.mtz", ldtype=2))
    with pytest.raises(ValueError, match="M/ISYM is 0"):
        validate_unmerged_mtz(_exported_mtz(tmp_path / "misym.mtz", misym=0))


def test_wavelength_band_plot(tmp_path):
    plot = wavelength_band_plot(_exported_mtz(tmp_path / "ok.mtz"))
    (data,) = plot["data"]
    assert data["x"] == [1, 101]
    assert data["y"] == pytest.approx([2.5, 3.5])
    assert data["error_y"]["array"] == pytest.approx([0.5, 0.0])
