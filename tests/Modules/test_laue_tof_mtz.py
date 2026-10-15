from __future__ import annotations

import gemmi
import numpy
import pytest

from xia2.Modules.Laue_TOF.laue_tof_mtz import (
    F_ALAMBD,
    F_DELAMB,
    F_PHIRANGE,
    I_LDTYPE,
    finish_unmerged_mtz,
    validate_unmerged_mtz,
)


def _exported_mtz(path):
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
                [1, 0, 0, 1, 1, 2e-3, 1e-4, 2.0],
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
        batch.ints[I_LDTYPE] = 2
        batch.floats[43] = 1.0
        for i, value in enumerate(umat.T.reshape(9)):
            batch.floats[6 + i] = value
        batch.dataset_id = 1
        mtz.batches.append(batch)
    mtz.write_to_file(str(path))
    return path


def test_finish_unmerged_mtz(tmp_path):
    mtz_file = _exported_mtz(tmp_path / "unmerged.mtz")
    with pytest.raises(ValueError, match="LDTYPE 3"):
        validate_unmerged_mtz(mtz_file, expected_batches=2)

    finish_unmerged_mtz(mtz_file, expected_batches=2)

    mtz = gemmi.read_mtz_file(str(mtz_file))
    assert [batch.ints[I_LDTYPE] for batch in mtz.batches] == [3, 3]
    assert [batch.floats[F_PHIRANGE] for batch in mtz.batches] == [1.0, 1.0]
    assert mtz.batches[0].floats[F_ALAMBD] == pytest.approx(2.5)
    assert mtz.batches[0].floats[F_DELAMB] == pytest.approx(0.5)
    assert mtz.batches[1].floats[F_ALAMBD] == pytest.approx(3.5)
    assert numpy.asarray(mtz.column_with_label("I").array) == pytest.approx(
        [2e3, 4e3, 3e3], rel=1e-5
    )
    assert numpy.asarray(mtz.column_with_label("SIGI").array) == pytest.approx(
        [1e2, 2e2, 1e2], rel=1e-5
    )


def test_validate_rejects_zero_misym(tmp_path):
    mtz_file = _exported_mtz(tmp_path / "unmerged.mtz")
    finish_unmerged_mtz(mtz_file, expected_batches=2)
    mtz = gemmi.read_mtz_file(str(mtz_file))
    data = numpy.array(mtz, copy=True)
    data[0, 3] = 0
    mtz.set_data(data)
    mtz.write_to_file(str(mtz_file))
    with pytest.raises(ValueError, match="M/ISYM is 0"):
        validate_unmerged_mtz(mtz_file)
