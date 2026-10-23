"""
Wavelength normalisation, scaling and merging with LAWLESS.

dials.scale and dials.merge are monochromatic: they have no per-observation
wavelength and no incident-spectrum model, so a Laue dataset cannot be scaled
with them. lawless, the Laue branch of aimless, fits the residual wavelength
dependence, scales and merges, reading dials.export's MTZ directly.

lawless is not in CCP4 yet, so it is looked up and reported rather than assumed,
and a run that cannot find it stops cleanly after the unmerged MTZ instead of
failing at the end of a long job.
"""

from __future__ import annotations

import logging
import os
import pathlib
import re
import shutil
import subprocess
import xml.etree.ElementTree as ElementTree

import gemmi
import numpy

from xia2.Driver.timing import record_step
from xia2.Handlers.Streams import banner
from xia2.Modules.Laue_TOF.laue_tof import LawlessParams
from xia2.Modules.Laue_TOF.laue_tof_mtz import check_wavelength_column
from xia2.Modules.Laue_TOF.util import result_files

xia2_logger = logging.getLogger(__name__)


def _candidates(names: tuple[str, ...]) -> list[pathlib.Path]:
    """Where a CCP4-side program might be, $PATH before the CCP4 installation."""
    found: list[pathlib.Path] = []
    for name in names:
        on_path = shutil.which(name)
        if on_path:
            found.append(pathlib.Path(on_path))
        ccp4 = os.environ.get("CCP4")
        if ccp4:
            candidate = pathlib.Path(ccp4) / "bin" / name
            if candidate.is_file():
                found.append(candidate)
    # Keep the order, drop repeats and links to the same file
    unique: list[pathlib.Path] = []
    for path in found:
        resolved = path.resolve()
        if resolved not in [other.resolve() for other in unique]:
            unique.append(path)
    return unique


def _ccp4_version(program: pathlib.Path) -> tuple[int, ...] | None:
    """The version a CCP4 program reports in its banner, or None."""
    try:
        result = subprocess.run(
            [os.fspath(program)],
            input="END\n",
            capture_output=True,
            encoding="utf-8",
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(
        r"version\s+(\d+)\.(\d+)\.(\d+)", result.stdout + result.stderr, re.I
    )
    if not match:
        return None
    return tuple(int(part) for part in match.groups())


def _find_program(
    names: tuple[str, ...], override: pathlib.Path | None, what: str
) -> pathlib.Path:
    """Where a CCP4-side program is, or a message saying how to say where."""
    if override and not override.is_file():
        raise FileNotFoundError(f"{what} is not at {override}")
    candidates = [override] if override else _candidates(names)
    for candidate in candidates:
        version = _ccp4_version(candidate)
        xia2_logger.info(
            f"Using {what} {_version_string(version)} from {candidate}"
            if version
            else f"Using {what} from {candidate}"
        )
        return candidate
    raise FileNotFoundError(
        f"Unable to find {what} ({' or '.join(names)}) on $PATH or in $CCP4/bin."
        f" Install it, or give its location with scaling.{names[0]}.executable."
    )


def _version_string(version: tuple[int, ...] | None) -> str:
    return ".".join(str(part) for part in version) if version else "of unknown version"


def _run(
    program: pathlib.Path,
    arguments: list[str],
    keywords: list[str],
    working_directory: pathlib.Path,
    logfile: pathlib.Path,
    name: str,
) -> None:
    """
    Run one of the CCP4 programs, keywords on stdin, output to a log file.

    They need the CCP4 environment for their libraries and symmetry data, so
    whatever is in the environment is passed through unchanged.
    """
    command = [os.fspath(program)] + arguments
    xia2_logger.debug(
        "Running %s with keywords:\n%s", " ".join(command), "\n".join(keywords)
    )
    with record_step(name):
        result = subprocess.run(
            command,
            cwd=working_directory,
            input="\n".join(keywords) + "\n",
            capture_output=True,
            encoding="utf-8",
        )
    logfile.write_text(result.stdout + result.stderr)
    result_files.record_log(logfile)
    if "Cannot open file" in result.stdout:
        raise ValueError(
            f"{name} could not open one of its input files. The command was:\n"
            + " ".join(command)
        )
    if result.returncode:
        raise ValueError(
            f"{name} returned error status {result.returncode}. Its log is"
            f" {logfile}, and it ended with:\n"
            + "\n".join(result.stdout.splitlines()[-15:])
        )


def _batch_range(mtz_file: pathlib.Path) -> tuple[int, int]:
    """The first and last batch number in an MTZ."""
    mtz = gemmi.read_mtz_file(os.fspath(mtz_file))
    numbers = [batch.number for batch in mtz.batches]
    if not numbers:
        raise ValueError(f"{mtz_file} has no batch records")
    return min(numbers), max(numbers)


def _wavelength_range(mtz_file: pathlib.Path) -> tuple[float, float]:
    """The range of the per-observation wavelength column in an MTZ."""
    label = check_wavelength_column(mtz_file)
    if not label:
        raise ValueError(f"{mtz_file} has no wavelength column")
    mtz = gemmi.read_mtz_file(os.fspath(mtz_file))
    values = numpy.asarray(mtz.column_with_label(label).array)
    return float(values.min()), float(values.max())


def lawless_keywords(
    hklin: pathlib.Path, params: LawlessParams, lorentz_applied: bool
) -> list[str]:
    """
    The keyword lines for one lawless run.

    The wavelength range and the batch range come from the data unless they were
    given.
    """
    first, last = _batch_range(hklin)
    lam_min, lam_max = _wavelength_range(hklin)
    if params.lam_min is not None:
        lam_min = params.lam_min
    if params.lam_max is not None:
        lam_max = params.lam_max
    lam_ref = (
        params.lam_ref if params.lam_ref is not None else 0.5 * (lam_min + lam_max)
    )

    keywords = [
        "TITLE xia2.laue_tof",
        f"RUN 1 BATCH {first} TO {last}",
        f"PROBE {params.probe.upper()}",
        f"SCALES {params.scales}",
    ]
    if params.normalisation == "chebyshev":
        keywords.append(
            f"LAUE NORMCHEBYSHEV {params.chebyshev_degree} {lam_min:.3f} {lam_max:.3f}"
        )
    elif params.normalisation == "gpr":
        keywords.append(f"LAUE NORMGPR {lam_min:.3f} {lam_max:.3f}")
        if params.gpr_bins:
            keywords.append(f"LAUE NORMGPRBINS {params.gpr_bins}")
    if params.normalisation != "none":
        keywords.append(f"LAUE NORMLAMREF {lam_ref:.3f}")
    if not lorentz_applied:
        xia2_logger.warning(
            "These data do not carry the Lorentz correction and lawless cannot"
            " apply it. The wavelength part of the factor will be absorbed by"
            " the normalisation curve, but sin^2(theta) will not, so the merged"
            " intensities will be mis-scaled between resolution shells."
            " Reintegrate with integration.lorentz=True."
        )
    if params.lambda_only:
        keywords.append("LAUE LAMBDAONLY")
    keywords.append("ANOMALOUS ON" if params.anomalous else "ANOMALOUS OFF")
    if params.sdcorrection:
        keywords.append(f"SDCORRECTION {params.sdcorrection}")
    if params.resolution:
        low, high = params.resolution
        keywords.append(f"RESOLUTION {low} {high}")
    keywords.append("OUTPUT MERGED UNMERGED")
    keywords += params.keywords
    return keywords


def run_lawless(
    working_directory: pathlib.Path,
    hklin: pathlib.Path,
    params: LawlessParams,
    lorentz_applied: bool,
    hklout_name: str = "scaled.mtz",
) -> pathlib.Path:
    """Normalise for wavelength, scale and merge."""
    program = _find_program(("lawless", "aimless"), params.executable, "lawless")
    working_directory = working_directory.resolve()
    # lawless writes LAMBDANORM, SCALES and ROGUES where it runs, and resolves
    # a relative HKLIN against somewhere else
    hklin = hklin.resolve()
    hklout = working_directory / hklout_name
    xia2_logger.notice(banner("Scaling with lawless"))  # type: ignore

    keywords = lawless_keywords(hklin, params, lorentz_applied)
    for line in keywords:
        xia2_logger.info(f"  {line}")

    _run(
        program,
        ["HKLIN", os.fspath(hklin), "HKLOUT", hklout_name, "XMLOUT", "lawless.xml"],
        keywords,
        working_directory,
        working_directory / "lawless.log",
        "lawless",
    )
    if not hklout.is_file():
        raise ValueError(f"lawless did not write {hklout}")

    result_files.record_data(hklout)
    result_files.record_data(working_directory / "lawless.xml")
    # The unmerged output carries the scaled intensities, with SCALEUSED as the
    # factor already applied, so nothing downstream may apply it again.
    unmerged = working_directory / f"{hklout.stem}_unmerged.mtz"
    if unmerged.is_file():
        result_files.record_data(unmerged)
    for extra in (
        "LAMBDANORM",
        "SCALES",
        "ROGUES",
        "ROGUEPLOT",
        "NORMPLOT",
        "CORRELPLOT",
        "ANOMPLOT",
    ):
        for path in sorted(working_directory.glob(f"{extra}*")):
            result_files.record_data(path)

    _report_lawless(
        working_directory / "lawless.xml", working_directory / "lawless.log"
    )
    return hklout


def lawless_xml(xmlout: pathlib.Path) -> ElementTree.Element | None:
    """
    The root of a lawless XML file, repaired if need be.

        lawless 0.0.1 opens <LAWLESS> and closes </AIMLESS>, so a strict parser
        refuses every file it writes.
    """
    text = xmlout.read_text()
    try:
        return ElementTree.fromstring(text)
    except ElementTree.ParseError:
        pass
    repaired = text.replace("</AIMLESS>", "</LAWLESS>")
    try:
        root = ElementTree.fromstring(repaired)
    except ElementTree.ParseError as e:
        xia2_logger.warning(f"Unable to read {xmlout.name}: {e}")
        return None
    xia2_logger.debug(
        f"{xmlout.name} closes <LAWLESS> with </AIMLESS>; read it as LAWLESS"
    )
    return root


NORMALISATION_DETAIL = (
    ("ReferenceWavelength", "reference wavelength"),
    ("LambdaMin", "from"),
    ("LambdaMax", "to"),
    ("Kernel", "kernel"),
    ("LengthScale", "length scale"),
    ("TrainingBins", "training bins"),
)


def lawless_log_sections(logfile: pathlib.Path) -> tuple[list[str], list[str]]:
    """The ASCII w(lambda) plot and the closing summary from a lawless log."""
    lines = logfile.read_text().splitlines()
    spectrum: list[str] = []
    summary: list[str] = []
    for i, line in enumerate(lines):
        if line.strip() == "w(lambda):":
            end = _first_line(lines, i, "marks reference wavelength")
            if end is not None:
                spectrum = lines[i : end + 1]
        elif "<!--SUMMARY_BEGIN-->" in line:
            end = _first_line(lines, i, "<!--SUMMARY_END-->")
            if end is not None:
                summary = lines[i + 1 : end]
    return spectrum, summary


def _first_line(lines: list[str], start: int, text: str) -> int | None:
    return next((j for j in range(start, len(lines)) if text in lines[j]), None)


def _report_lawless(xmlout: pathlib.Path, logfile: pathlib.Path) -> None:
    """Log lawless's warnings, its w(lambda) plot and its summary."""
    _report_lawless_xml(xmlout)
    if not logfile.is_file():
        return
    for section in lawless_log_sections(logfile):
        if section:
            xia2_logger.info("\n".join(section))


def _report_lawless_xml(xmlout: pathlib.Path) -> None:
    if not xmlout.is_file():
        xia2_logger.warning(f"lawless wrote no {xmlout.name}")
        return
    root = lawless_xml(xmlout)
    if root is None:
        return

    for warning in root.iter("WarningMessage"):
        for line in (warning.text or "").splitlines():
            if line.strip():
                xia2_logger.warning(f"lawless: {line.strip()}")

    for tag in ("WavelengthNormalisationGPR", "WavelengthNormalisationChebyshev"):
        block = root.find(tag)
        if block is None:
            continue
        detail = [
            f"{label} {block.find(name).text.strip()}"  # type: ignore[union-attr]
            for name, label in NORMALISATION_DETAIL
            if block.find(name) is not None and block.find(name).text  # type: ignore[union-attr]
        ]
        if detail:
            xia2_logger.info(
                f"Wavelength normalisation ({tag.replace('WavelengthNormalisation', '')}):"
                f" {', '.join(detail)}"
            )
        break


def _report(
    working_directory: pathlib.Path,
    scaled: pathlib.Path,
    unmerged_mtz: pathlib.Path,
) -> None:
    """
    Write the run report from the scaled data.

        A failure here is reported but does not fail the run: the data it describes
        is already written.
    """
    from xia2.Modules.Laue_TOF.laue_tof_report import generate_report

    # lawless writes the unmerged companion beside the merged output
    scaled_unmerged = scaled.with_name(f"{scaled.stem}_unmerged.mtz")
    if not scaled_unmerged.is_file():
        xia2_logger.warning(
            f"lawless wrote no {scaled_unmerged.name}, so there is nothing"
            " unmerged to report on. Add OUTPUT UNMERGED to"
            " scaling.lawless.keywords."
        )
        return
    try:
        generate_report(
            working_directory,
            scaled_unmerged,
            unmerged_mtz=unmerged_mtz,
            lawless_directory=scaled.parent,
        )
    except Exception as e:
        xia2_logger.warning(f"Unable to write the report: {e}")


def scale(
    working_directory: pathlib.Path,
    unmerged_mtz: pathlib.Path,
    lawless_params: LawlessParams,
    lorentz_applied: bool,
    steps: list[str],
) -> pathlib.Path | None:
    """
    Run the scaling tail, as far as the steps and the programs available allow.

    A missing program is reported and the run stops here with the unmerged MTZ
    in hand, rather than failing after everything that came before it.
    """
    if "lawless" not in steps:
        return None
    try:
        scaled = run_lawless(
            working_directory,
            unmerged_mtz,
            lawless_params,
            lorentz_applied,
        )
        if "report" in steps:
            _report(working_directory, scaled, unmerged_mtz)
        return scaled
    except FileNotFoundError as e:
        xia2_logger.warning(
            f"{e} The unmerged data is in {unmerged_mtz}, so the scaling can be"
            " run separately once the program is available."
        )
        return None
