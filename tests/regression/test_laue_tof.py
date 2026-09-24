from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import h5py
import numpy
import pytest
from dxtbx.serialize import load

from xia2.Modules.Laue_TOF.laue_tof import (
    FileInput,
    check_inputs_readable,
    readable_images,
)

# The MANDI test data is a 3 GB event file, too big for dials_data, so the
# MANDI tests run only where a copy exists. Default location is
# work/mandi/DataFiles in the repository (gitignored, since /work is);
# XIA2_LAUE_TOF_MANDI overrides it.
MANDI_DEFAULT = (
    pathlib.Path(__file__).parents[2]
    / "work"
    / "mandi"
    / "DataFiles"
    / "MANDI_13378.nxs.h5"
)


@pytest.fixture
def mandi_binned() -> pathlib.Path:
    """A MANDI file holding TOF bins, or skip if there is no copy of it here."""
    path = pathlib.Path(os.environ.get("XIA2_LAUE_TOF_MANDI", MANDI_DEFAULT))
    if not path.is_file():
        pytest.skip(
            f"MANDI test data not found at {path}. Copy MANDI_13378.nxs.h5"
            " there, or point XIA2_LAUE_TOF_MANDI at it, to run the MANDI"
            " tests."
        )
    if not readable_images(path):
        pytest.skip(
            f"{path} holds raw event data, which this pipeline does not reduce."
            " Bin it to time-of-flight first."
        )
    return path


def test_unreduced_data_is_not_readable(tmp_path):
    """
    A file a format class claims but cannot read an image from is not usable.

    dxtbx's understand() answers which instrument wrote a file, not whether the
    data in it have been reduced: FormatMANDI claims anything whose instrument
    name is MANDI, and FormatESSNMX anything named NMX, raw or not. So what
    decides whether this pipeline can process a file is whether an image can
    actually be read from it.
    """
    mandi = tmp_path / "raw_mandi.h5"
    with h5py.File(mandi, "w") as fh:
        entry = fh.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry.create_group("instrument").create_dataset(
            "name", data=numpy.array([b"MANDI"])
        )
    assert readable_images(mandi) == 0

    nmx = tmp_path / "raw_nmx.h5"
    with h5py.File(nmx, "w") as fh:
        entry = fh.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        entry.create_group("instrument").create_dataset(
            "name", data=numpy.bytes_("NMX")
        )
    assert readable_images(nmx) == 0

    not_data = tmp_path / "notdata.txt"
    not_data.write_text("not a data file")
    assert readable_images(not_data) == 0

    with pytest.raises(ValueError, match="cannot read"):
        check_inputs_readable(FileInput(images=[os.fspath(mandi)]))


def test_binned_data_is_readable(mandi_binned):
    """A file holding TOF bins reports its images and passes the check."""
    assert readable_images(mandi_binned) > 1
    check_inputs_readable(FileInput(images=[os.fspath(mandi_binned)]))


def test_import(regression_test, mandi_binned, tmp_path):
    """Import a binned MANDI run, which is what the pipeline starts from."""
    result = subprocess.run(
        [
            shutil.which("xia2.laue_tof"),
            f"image={os.fspath(mandi_binned)}",
            "workflow.steps=find_spots",
        ],
        cwd=tmp_path,
        capture_output=True,
        encoding="utf-8",
    )

    imported = tmp_path / "orientation_1" / "imported.expt"
    assert imported.is_file(), result.stdout
    expts = load.experiment_list(imported, check_format=False)
    assert expts.all_tof()
    assert len(expts[0].detector) > 1
    assert expts[0].scan.get_num_images() == readable_images(mandi_binned)
