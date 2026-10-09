# xia2.laue_tof — neutron Laue time-of-flight pipeline

## Status (2026-10-09)

Runs end to end: import → find_spots → index → refine → tof_integrate → combine
→ export → unmerged MTZ → pointless → lawless → report. What is left is tuning,
tests and docs, not structure.

| file | what it is |
|---|---|
| `Modules/Laue_TOF/laue_tof.py` | `*Params` dataclasses, `setup_from_phil`, `run_xia2_laue_tof` |
| `Modules/Laue_TOF/laue_tof_phil.py` | the phil scope |
| `Modules/Laue_TOF/data_integration.py` | the step functions and `run_data_integration` |
| `Modules/Laue_TOF/laue_tof_mtz.py` | unmerged Laue MTZ writer + `validate_unmerged_mtz` |
| `Modules/Laue_TOF/laue_tof_scale.py` | pointless and lawless steps |
| `Modules/Laue_TOF/laue_tof_report.py` | the report step |
| `Modules/Laue_TOF/util.py` | `report_timing` (copied, never imported from SSX) |
| `cli/laue_tof.py` | the command, registered in `setup.py` |

### Done

| commit | date | what |
|---|---|---|
| `efd855e9` | 10-09 | `d_min`/`d_max` redefined as a cut on integrated reflections |
| `48e54cc5` | 10-08 | whole run processed together; `joint_indexing=False`; refine at `hierarchical`/level 1 |
| `816b6d56` | 10-08 | `max_cell` always fixed; `max_cell_fallback` = 150 Å |
| `000554d0` | 09-26 | ISAW `.DetCal` calibration at import |
| `d62b2300` | 09-24 | integrate to the geometric limit; Lorentz hard-coded on |
| `94ad4a4b` `863a78db` `1e7cce04` | 09-21 | scaling tail, unmerged MTZ, report |

### Next, in order

1. **Rank indexing solutions on more than RMSD.** `data_integration.py` keeps the
   best accepted attempt by median RMSD; the study calls for the geometric mean of
   indexed count, fraction and RMSD. The 50-bin rubredoxin case is why: `fft1d`
   finds 1.57 px on ~1300 of 9102 spots and `min_indexed_fraction` (0.25) rejects
   it, the denominator being mostly junk. Probably means judging the fraction over
   a resolution range the data support.
2. **Extend the ladder past `indexing.method`** to the threshold algorithm and
   `min_spot_size`. A `min_spot_size` rung needs `dials.find_spots` re-run, which
   the ladder does not currently do.
3. **Warn when a cell is up against `max_cell`.** Fixing the default moved the
   failure rather than removing it: a cell longer than the fallback still
   truncates silently.
4. **Tests, then docs.** 3 tests, all gated on a MANDI file that may be absent; no
   SXD case, no converter unit test, no `doc/sphinx` page. The multi-exposure
   `imageset_id` bookkeeping is tested only by hand.
5. **Re-add TOF binning** as a step that does not classify input on the pipeline's
   behalf.

### Testing convention

Pass the cell and symmetry; it is rare not to know them in a Laue-TOF experiment.
For the NMX rubredoxin simulations the truth is in the McStas input
`Rubredoxin_solvent.lau`: `P 21 21 21`, `34.375 35.183 43.994 90 90 90`. PPase is
`C 1 2 1 106.1 95.5 113.7 90 98.1 90`; MANDI CuZnSOD is
`P 65 2 2 65.73 65.73 151.88`. It matters: blind, the 50-bin rubredoxin data
returns a halved cell (a = 17.2 Å) at 8.42 px; with the cell given, the full
`34.70 35.56 44.42` at 1.57 px.

## Context

xia2 had no neutron Laue-TOF support: existing pipelines assume one scalar
wavelength per sweep and `load_imagesets` keeps only rotation sequences. Laue-TOF
is polychromatic — each reflection's wavelength comes from its arrival time, the
beam is a `PolychromaticBeam`, and intensities need normalising by an incident
spectrum (vanadium), background (empty run) and absorption.

**Experiment shape**: one crystal, several **stationary** exposures at different
settings. Each exposure is a full Laue-TOF measurement (a stack of TOF bins, one
`ImageSequence`) with the crystal still — the contrast with rotation
crystallography. The whole run goes through the DIALS steps together, then combine
→ unmerged MTZ → POINTLESS → LAWLESS. Instruments: ISIS SXD (the only public
`dials_data`), MaNDi, ESS NMX.

DIALS supplies the TOF toolchain; this wires it into one command. Scaling is done
by two locally built CCP4-side programs, POINTLESS 1.16.1 and LAWLESS 0.0.1 (the
Laue branch of aimless).

## Decisions

- No `dials.scale`/`dials.merge` — monochromatic only. Lawless replaces them.
- Vanadium and empty runs optional but recommended; corrections skipped if absent.
- SHELX HKLF2 of the unscaled combined data is a checkpoint; lawless's scaled MTZ
  is the deliverable.
- Indexing defaults come from the NMX study, not DIALS' own, and indexing is a
  **retry ladder with an acceptance test**, not one call.
- Formats: `FormatISISSXD`, `FormatMANDI`, `FormatESSNMX`. The `*Laue` variants are
  for TOF-flattened data and are out of scope.
- `refine_bravais_settings` stays a toggleable, informational step.
- Local pointless/lawless builds, not CCP4 9's — paths are phil-configurable.

## DIALS facts worth knowing

- `ExperimentList.all_tof()` / `all_laue()` identify this data. Assert `all_tof()`
  after import.
- **Avoid `FormatMANDILaue` / `FormatESSNMXLaue`**: they read the `laue_data` group,
  flattened in TOF, which is a different problem. Their `understand()` keys on
  `"laue_data" in …detector_panel_0`, so a file holding both groups can be claimed
  by the wrong one. `FormatESSNMX_McStas` is for the superseded McStas format;
  current McStas output is identical to real NMX data.
- `dials.find_spots` auto-selects `TOFSpotFinder`.
- `dials.tof_integrate` options live under `corrections{}`, not top level, and
  `corrections.lorentz` **defaults to False**.
- **`dials.export format=mtz` is refused for TOF data**, and would be wrong anyway:
  no per-reflection wavelength, `M/ISYM = 0`, `LDTYPE = 2`. Hence xia2's own writer.
- `dials.export format=shelx` takes the batch from `imageset_id`, which is 0 in
  every per-exposure file — the combine step renumbers it.
- `dials.tof_integrate` leaves **no record of the Lorentz correction**: it
  multiplies intensities in C++ and writes no column or flag, only a log line. So
  the integrate step stamps a `lorentz_applied` column, which survives combine and
  export.

## Gotchas found on real data

Each cost a session; none is visible from the code.

- **`joint_indexing` must be passed explicitly.** A TOF experiment has a goniometer
  and a scan over its bins, so `is_still()` is False; with >1 experiment
  `joint_indexing=Auto` takes the all-rotation branch and aborts. Hard-coded False.
- **`dials.index` output cannot be read positionally.** With `joint_indexing=False`
  it indexes in a process pool and extends the experiment list inside an
  `as_completed` loop, so order is completion order, and failures are dropped — the
  list is shorter than the input, in an arbitrary order. `imageset_id` is reliable
  (DIALS resets it to the original index); everything maps back through it.
- **A rung that solves nothing writes a column-less reflection table**, which raises
  `DXTBX_ASSERT(it != end())` when asked for a flag.
- **`dials.tof_integrate` assumes one detector shape and TOF binning** for the set,
  reading `image_size`/`tof_size` from `experiments[0]`. Two reductions of one run
  at 50 and 100 bins import together without complaint, so import rejects that.
- **`panels=hierarchical` is a no-op without a level.** The NMX hierarchy is flat
  (three panels on one root group), so level 0 is the whole detector as one body —
  16 parameters, identical residuals to `automatic`. Level 1 reaches the panels, 28
  parameters, per-panel RMSD 0.32 px.
- **`max_strong` is per exposure**, since the predictor assertion it guards is
  reached per refinement job (`dials.refine` refines each exposure separately and
  recombines).
- **`spotfinding.d_min`/`d_max` cannot exist.** A DIALS resolution cut is a pixel
  mask built from a `ResolutionMaskGenerator`, which needs a fixed `s0` that a
  `PolychromaticBeam` lacks (`PolychromaticBeam has no fixed s0`, `beam.h(741)`),
  and a Laue-TOF pixel's resolution depends on the wavelength that arrived there.
  Top-level `d_min`/`d_max` therefore cut the **integrated reflections** in
  `combine`, by their own `d`, before export, MTZ and scaling. Deleting the names
  instead was worse: phil resolves an unambiguous suffix, so `d_min=1.8` silently
  became `integration.calculated_d_min=1.8`.
- **`max_cell` must never be left to DIALS.** On rubredoxin the nearest-neighbour
  estimate gave 23.3 Å for a 44.1 Å axis and returned a halved cell that passed
  every acceptance check — right space group, 0.59 px, 31 % indexed. It is always
  fixed now: `max_cell_multiplier` × the longest edge of a given cell, else
  `max_cell_fallback` = 150 Å. Verified at 150 on rubredoxin (34×35×44) and PPase
  (primitive 71.7×71.7×114.5, where the estimate happens to land at 149.4 Å).
- **MANDI's `c` is 151.88 Å**, longer than the fallback — it supplies a cell, so the
  derived rule covers it, but a MANDI-like run indexed blind would not be.
- **`FormatMANDI` hard-codes a 2.0–4.0 Å band** where the histogram gives 1.92–4.08,
  which alone moves its geometric `d_min` from 1.066 to 1.027 Å. So the wavelength
  range is derived from the scan, not from `beam.get_wavelength_range()`.
- **`min_spot_size` 15 is tuned for NMX and is too big for MANDI**, where it
  discards 59 % of 558 spots. On NMX it discards 99 % of 187 300 and leaves 2296,
  which is not a problem. The step reports the fraction and advises only when few
  spots survive.
- **TOF binning changes what indexes.** Blind, rubredoxin gives the full cell at 100
  bins and a halved one at 50 and 200 — `max_cell` does not rescue it, a known cell
  does. Binning is a parameter, not a default.
- **Identical `image=` values are deduplicated by phil**, and `resolve_paths`
  resolves symlinks, so two links to one file become one path but still two entries.
  Use real copies to fake a multi-exposure test.

### Two DIALS defects

1. **`known_symmetry.A_matrix` corrupts a centred cell.** The indexer expects a
   primitive-setting crystal; the conversion in `known_orientation.py` assigns to
   the loop variable and is discarded, while `_apply_symmetry_post_indexing`
   converts an already-centred cell again (`C 1 2 1 106.1 96.1 114.4 · 98.1` →
   `184.1 158.6 102.1 · 90`). Convert in the caller, which makes the missing
   conversion a no-op; the acceptance test's cell check would catch a fixed DIALS
   converting twice.
   ```python
   cb_op = crystal.get_space_group().info().change_of_basis_op_to_primitive_setting()
   A = crystal.change_basis(cb_op).get_A()
   ```
2. **TOF predictor assert on large spot lists** — around 40 000 reflections trips
   `DIALS_ASSERT(table.nrows() == h.size())` in `reflection_predictor.h` during
   refinement. Hence `max_strong`.

Also unfixed: the TOF profile fit runs in a `ThreadPool` worker that catches
nothing, so any throw kills the process rather than one reflection. That is what
the ellipse-mask defect (`compute_weighted_ellipsoid`, "Values must not all be
zero") rides on — use `mask=seed_skewness` on sparse data.

## Indexing defaults (from the NMX study)

Source: `~/ppase_allconfigs/NMX_indexing_{recommendations,study}.html` (PPase,
`C 1 2 1`, 4172 trials over eleven simulated NMX layouts at two binnings, DIALS
3.30 + optuna), plus `~/ppase_allconfigs/cco` (CcO, `P 21 21 21`, 182 204 178).
Per-run recipes in `configN/opt_NNbins/best*/`; every trial in `results.jsonl`.
The search took datasets that autoindex from 4 of 22 to 15 of 22, and best RMSD
from 2.3 px to 0.59 px. **These numbers were measured on 3.30, not dials-dev.**

| parameter | DIALS default | here | why |
|---|---|---|---|
| `min_spot_size` | 6 | **15** | the largest single improvement: RMSD 2.3 → 0.6–0.8 px, at some cost in count |
| `threshold.algorithm` | `radial_profile` | **`dispersion_extended`** | wins on four of six near PPase layouts; `radial_profile` wins more often on CcO, so a default, not a rule |
| `indexing.method` | `fft3d` | **ladder** `fft3d` → `real_space_grid_search` → `fft1d` | `fft1d` was the finding: the only method that indexed the 1 m layouts, its 1-D projection not needing a populated 3-D grid |
| `max_cell` | Auto | **always fixed** (see Gotchas) | NN estimation is unreliable on sparse spot lists; every difficult-layout success used a fixed value |
| `outlier.algorithm` | `auto` | **`mcd`**, iqr 0.8, `separate_panels=True` | distant layouts carry a heavy residual tail (p99/p50 47.6 against 8.9) |
| `d_min_start` | Auto | **4.0** | helps the difficult layouts, immaterial for the near ones |

Consequences:

- **A ladder, not a recipe.** Near the boundary, success is a rare event rather than
  a property of the data: config8 landed 2 acceptable solutions in 212 trials,
  config11 4 in 190. One recipe will fail on indexable data.
- **Accept on more than "it indexed".** Maximising indexed fraction alone picks a
  small clean spot list (88.7 % of 909 beating 72.1 % of 2837), and an a↔b swap can
  index 1988 reflections on the wrong cell. Rank on count, fraction and RMSD
  together — **still outstanding**, item 1.
- **Seed from a sibling.** Every exposure is the same crystal, and six NMX
  configurations indexed independently agreed on orientation to 0.04–0.09°. Seeding
  also rescues layouts that never autoindex: config7 and config9 gave nothing in
  480 and 330 de-novo trials but process normally when seeded. So "did not index"
  is a routine outcome — seed, log it, carry on.
- **Binning is a parameter.** PPase indexes better coarse (50) and integrates
  better fine (200); large-cell CcO prefers 200 for both (indexed fraction
  0.70–0.84 against 0.24–0.58). A layout that fails at one binning is worth
  retrying at the other.

## Resolution: observed, then calculated

`integration_type=observed` integrates only the indexed spots — a few hundred of
the tens of thousands the geometry reaches. `calculated` integrates every predicted
reflection out to a limit, defaulting to the geometric one: `d = λ / (2 sin θ)` at
the shortest wavelength and largest scattering angle any pixel sees, which nothing
can be recorded beyond.

| | wavelengths | 2θ max | geometric d_min | predicted | observed |
|---|---|---|---|---|---|
| NMX config1, 50 bins | 1.87–3.54 Å | 134.0° | **1.016 Å** | 93 304 | 1359 |
| MANDI CuZnSOD, 50 bins | 1.92–4.08 Å | 139.3° | **1.027 Å** | 82 630 | ~230 |

- **Both passes always run** under `calculated`: the observed pass chooses the
  method and mask at a hundredth of the cost of comparing them over 93 000
  shoeboxes, and its output stays on disk as a checkpoint.
- **Cost guard `max_predicted`** (250 000), counted with `TOFReflectionPredictor`
  before any shoebox is extracted, and a run over it stops with the d_min that
  would bring it under. Found by prediction, not from `d^-3`: a Laue-TOF experiment
  does not record a sphere, so the analytic estimate is optimistic (2.69 Å against
  a true 3.31 Å).
- **The wavelength band bounds the calculated pass only.** An observed reflection is
  inside it by construction; a predicted one outside it would be integrated as
  noise. This is what protects MANDI, whose stated band is not its measured one.
- **Profile fitting a calculated pass is not worth it on MANDI**: `profile_1d_ibix`
  fits 98.9–100 % of observed spots but 55–61 % of everything predicted, so the file
  mixes two populations (42 887 observations with no profile value) and summation
  wins on every statistic. `min_profile_fraction` does not catch it — 58 % passes.
  Use `method=summation` for calculated passes until the 1D profile parameters are
  tuned.
- Cut afterwards on CC½ in lawless. Re-merging costs nothing; re-integrating costs
  everything.

## Integration choices

- `dials.tof_integrate` always writes `intensity.sum.*`; a profile method adds
  `intensity.prf.*` beside it. So "fall back to summation" is only for an outright
  failure, and both intensities reach the MTZ (`I`/`SIGI` and `IPR`/`SIGIPR`), a
  failed fit getting the missing-number flag, never a zero.
- `profile_1d_ibix` used to abort the run on small intensities
  (`DIALS_ASSERT(A >= min_bounds[0] …)`: `A` is fixed at 1.0 in normalised units
  while `max_bounds[0]` is `1e4 · intensity_max` in raw units, so the data's *scale*
  decided). **Fixed by [#3281](https://github.com/dials/dials/pull/3281)**, which
  bounds `A` by phil and corrects an inverted FWHM guard — that second half lifted
  the fit fraction on config3 from 11 % to 62 %. The fix is in `dials-dev` only.
- A profile run finishing is not the method working: `profile_3d_gutmann` completed
  having fitted 13 of 1359 reflections. Hence `min_profile_fraction` (0.5), below
  which summation is used.
- **The two masks give different data, not different metrics of it** (summation, one
  NMX orientation): `ellipse` 7.0 s wall, mean I/σ 36.5, 71.6 foreground px;
  `seed_skewness` 26.5 s, 59.6, 322.2 px. CC between the intensities 0.95. Higher
  I/σ is not proof of better, so `mask=both` integrates each way and keeps the
  winner by `mask_metric`, and the run is integrated in one call so every batch is
  treated alike.
- `dials.export format=shelx` refuses its own `intensity=auto` when both columns are
  present, so the export step always passes an explicit `intensity=`.

## Lorentz: exactly once

`corrections.lorentz=True` is hard-coded and the phil option is gone. The lawless
branch in use has no `LORENTZ` keyword and aborts if given one, so integration is
the only place the factor can be taken.

Measured on MANDI CuZnSOD. `L = sin²θ/λ⁴`, and the halves differ:

- `λ⁴` is **degenerate with the wavelength normalisation** — lawless says so in
  `mtz_unmerge_io.cpp` — so it makes no difference whether it is applied or absorbed
  by `w(λ)`. The MaNDi reduction leaves it in deliberately, which is why its
  normalisation curve rises where ours falls.
- `sin²θ` is **not** degenerate: it sets relative scale between resolution shells,
  and nothing downstream supplies it. Applying it wins 12 of 20 shells and the outer
  shell (CC½ 0.500 against 0.457). Overall CC½ looks worse (0.786 against 0.927)
  purely because the correction re-weights the shells.

Consequence: intensities land ~1e-3 of counts, which makes **pointless fail
outright** (`NormaliseInRuns::apply Unused run 1`). The MTZ writer therefore scales
`I`/`SIGI` and `IPR`/`SIGIPR` by a power of ten to a counts-like median and logs the
factor — a constant on both columns is absorbed by the scale factors and changes no
statistic.

## Scaling tail (POINTLESS → LAWLESS)

Reference: `~/aimless/lawless.md`, `~/aimless/CLAUDE.md`.

`dials.tof_integrate` normalises by the measured incident spectrum, which leaves
per-orientation scale, relative B, absorption and any **residual wavelength
dependence**. Lawless fits those: aimless plus a wavelength normalisation `w(λ)`,
as a Chebyshev polynomial or a Gaussian process, with `PROBE NEUTRON` defaults.
Pointless comes first because lawless needs reflections sorted and on a consistent
setting.

### The unmerged MTZ

Written with gemmi by `laue_tof_mtz.py`. Every requirement below has been observed
to fail **silently**, so `validate_unmerged_mtz` re-reads the file and fails loudly
on each.

| item | requirement | failure if wrong |
|---|---|---|
| `LAMBDA` | per-observation wavelength from `wavelength_cal` | lawless stops when `LAUE` is given |
| `M/ISYM` | present and **1 throughout** | pointless reads the file as merged, collapses observations, drops columns |
| batch records | one per exposure, `LDTYPE = 3`, a cell, non-zero `PHIRANGE`, `BSCALE = 1`, `UMAT` = refined orientation | no records ⇒ pointless reports one batch regardless of `BATCH`; no `UMAT` ⇒ `SECONDARY`/`ABSORPTION` refines flat with no warning |
| datasets | `HKL_base` as dataset 0, data in ≥ 1, batch `LBSETID` pointing at the data set | `Dataset::pxdname, setid not found` |
| `BATCH` | one per exposure, **1-based** (SHELX batches are 0-based) | mis-scaled or merged-across-settings data |
| `ROT` | per-batch φ, optional | lawless substitutes the batch number; harmless for stationary exposures |

Batch word offsets (`I_LDTYPE=14`, `F_UMAT=6`, `F_PHISTT=36`, `F_PHIRANGE=47`,
`F_BSCALE=43`, `F_ALAMBD=86`, `F_DELAMB=87`) are copied from
`~/LADI_files/TIM_2020/aaron/addbatches.py`, not re-derived.

### Pointless

```
pointless HKLIN unmerged.mtz HKLOUT sorted.mtz XMLOUT pointless.xml <<EOF
LAUEGROUP P 6/m m m
CHOOSE SPACEGROUP P 65 2 2
EOF
```

- **Always give the symmetry.** Pointless cannot determine the Laue group from
  unnormalised Laue data: the wavelength dependence inflates disagreement between
  mates, so even the identity scores badly. Treat any suggestion as informational;
  never let it silently reindex.
- **Chosen by version, not by path.** CCP4 9 ships 1.13.6, which drops the
  wavelength column; the local build is 1.16.1. Every candidate on `$PATH` and in
  `$CCP4/bin` is asked its version, the first newer than 1.13.6 wins, and the
  wavelength column is checked afterwards as a backstop.

### Lawless

```
lawless HKLIN sorted.mtz HKLOUT scaled.mtz XMLOUT lawless.xml <<EOF
PROBE NEUTRON
RUN 1 BATCH 1 TO <nbatches>
SCALES BATCH BFACTOR ON SECONDARY 0
LAUE NORMGPR <lam_min> <lam_max>
LAUE NORMLAMREF <lam_ref>
SDCORRECTION NOREFINE 1.0 0.0 0.0
REJECT 4
OUTPUT MERGED UNMERGED
EOF
```

- `PROBE NEUTRON` switches off automatic `ANOMALOUS`, zeroes polarisation, makes
  `SCALES BATCH` the default (φ is not a scaling variable here) and reports
  rejections against 2θ.
- `NORMGPR` is the better default: it follows structure a low-order polynomial
  cannot and is **fixed after its pre-pass**, so it cannot trade against the other
  scale terms. `NORMCHEBYSHEV` is refined jointly. They are **mutually exclusive —
  giving both is fatal**. Band limits come from the `wavelength_cal` range.
- **Both programs need absolute HKLIN paths.** Lawless writes its outputs into its
  working directory, so it gets one of its own, and a relative HKLIN then resolves
  against the wrong place — reported, misleadingly, as *"LAUE requires a wavelength
  for every observation"*. A `Cannot open file` is now caught directly.
- **Lawless 0.0.1 writes invalid XML**: it opens `<LAWLESS>` and closes
  `</AIMLESS>`, left from the rename, so strict parsers refuse it. The reader
  repairs that one mismatch. Worth fixing in lawless.
- Statistics live in `Result/Dataset/<name>/{Overall,Inner,Outer}`, and the fitted
  curve in `WavelengthNormalisationGPR` or `…Chebyshev`.
- Outputs registered: `scaled.mtz`, the unmerged output, `lawless.xml`, the log, and
  **`LAMBDANORM`** — a self-contained gnuplot script of `w(λ)` (GP carries a 1σ
  band, Chebyshev does not).
- **One batch per exposure** means few parameters per batch, so with a handful of
  settings `SCALES BATCH BFACTOR ON SECONDARY 0` may be over-parameterised. Keep
  `scaling.lawless.scales` overridable and report the parameter/observation ratio.

## The report

`laue_tof_report.py` produces `xia2.laue_tof-report.{html,json}` — xia2's own report
page with one panel added, not a new kind of report. Built from
`Report.from_unmerged_mtz` on lawless's `scaled_unmerged.mtz` and rendered in
process, so a Laue user gets the tables and plots they already know.

The added panel holds what the standard page cannot know: the fitted `w(λ)` with its
1σ band and the binned observations behind it, read from `LAMBDANORM`'s inline
blocks; each exposure's wavelength band from `ALAMBD`/`DELAMB`; profile against
summation intensities from the **unscaled** MTZ, the only file carrying both; and
what lawless fitted, from the XML.

`templates/report.html` gained an empty `{% block extra_panels %}` for this — one
line upstream, no duplicated template, reusable by any pipeline. A failure writing
the report does not fail the run.

## TOF binning (removed, to be re-added)

`dials.import` needs a TOF-binned file; raw event data must be reduced first, by
**essnmx** (`~/ess/packages/essnmx`). The in-pipeline `bin` step was removed on
2026-09-23 because it tangled binning with reading and the MANDI route mutates the
user's raw files. The pipeline now only checks that dxtbx can read each input and
names any file it cannot. Old code is in the history.

```
essnmx-reduce   --input-file raw.h5             --nbins 50 --output-file binned50.h5
essmandi-reduce --input-file mandi_event.nxs.h5 --nbins 50 --output-file binned50.h5
```

When re-added it must not classify input on the pipeline's behalf. Note
`ess/nmx/nexus.py::_set_default_instrument` writes `entry/instrument/name = "NMX"`
unconditionally, so a MANDI reduction is claimed by `FormatESSNMX`, whose
`get_detector` then raises `KeyError: 2` — fix that in essnmx first.

MANDI is histogrammed in place by dxtbx (bin edges and counts written into the event
file, events kept). Two dxtbx problems worth fixing:
`generate_histogram_data_for_panel` builds a pool task per pixel, ~1 min per panel
and ~40 min for a 40-panel run, where the same histogram is a few seconds of
`np.bincount`; and `get_time_range_for_panel` returns
`event_time_zero + event_time_offset`, which for a 23 h run gives 14700–113000 µs
where every event arrives between 14739 and 31407.

## Environment

```
source .venv/bin/activate
export PATH=~/dials/dials-dev/conda_base/bin:$PATH
```

The second line matters: `.venv/bin` holds the `xia2.*` scripts but no DIALS
dispatchers, and every step resolves its program with `shutil.which`.

| what | where |
|---|---|
| DIALS (in use) | `~/dials/dials-dev` — branch `tof-profile-1d-mantid`, `7e37c2844`. Has #3281, and adds `profile_1d_mantid` |
| DIALS 3.30 (**do not use**) | `~/dials/dials-v3-30-0` — clean at `main`, but its built extensions come from newer source, so `profile_1d_ibix` fails a Boost.Python signature check and integration silently falls back to summation. The build the NMX study used |
| essnmx / essreduce | `~/ess/packages/` (editable; scipp from PyPI, numpy pinned to the DIALS build's) |
| xia2 (this fork) | editable, so `xia2.*` runs this tree and not the DIALS tree's `modules/xia2` |
| pointless 1.16.1 | `~/pointless/pointless` |
| lawless 0.0.1 | `~/aimless/build/aimless` |
| NMX indexing study | `~/ppase_allconfigs/` (+ `cco/`) |

Two venv pieces are hand-made and **will not survive recreating it**:
`.venv/share/cctbx`, a symlink to the DIALS `conda_base/share/cctbx` (without it
`libtbx.load_env` raises `FileNotFoundError`, because libtbx looks under
`sys.prefix`, which in a venv is the venv), and `.venv/bin/pytest`, a shim running
this venv's interpreter. `ruff` is not in the venv;
`~/miniforge3/envs/nmx_strategy/bin/ruff` is the one used.

To rebuild against another DIALS tree: create the venv from that tree's
`conda_base/bin/python` with `--system-site-packages`, remake those two, then
`pip install -e .` here and `-e` the two `ess/packages` — all under a
`numpy==<that build's version>` constraint, or pip moves numpy and the compiled
extensions stop loading. An editable `dxtbx` is unnecessary when the tree's
`conda_base` already resolves `dxtbx` to its own `modules/dxtbx/src`.

## Verification

1. **Smoke**: `xia2.laue_tof -h`; pointless reports 1.16.1; lawless reports 0.0.1.
2. **End to end on SXD** (N=1, the only public data):
   ```
   dials.data get isis_sxd_example_data
   xia2.laue_tof image=<…>/sxd_nacl_run.nxs \
       unit_cell=5.64,5.64,5.64,90,90,90 space_group=225 \
       vanadium_run=<…>/sxd_vanadium_run.nxs empty_run=<…>/sxd_empty_run.nxs
   ```
   Expect the full chain to `scaled.mtz` + `LAMBDANORM`, `all_tof()` true on
   `imported.expt`, and `wavelength_cal` on the integrated reflections.
3. **MTZ**: `validate_unmerged_mtz` passes with one batch per exposure.
4. **Pointless passthrough**: `sorted.mtz` still has `LAMBDA` and the observation
   count matches — a large drop means it was read as merged.
5. **Lawless**: no "no wavelength column" stop; `PROBE NEUTRON` block printed;
   `LAMBDANORM` close to flat, vanadium normalisation having already been applied.
6. **Indexing defaults**: on `~/ppase_allconfigs/configN/`, the default recipe
   reproduces the study's balanced pick to within noise, and the ladder recovers a
   layout only `fft1d` indexes. A seeded C-centred crystal comes back with the input
   cell, not a doubled one.
7. **Multi-exposure**: several `image=`; one batch per exposure through to the MTZ,
   and `RUN 1 BATCH 1 TO N` covering them all.
8. **Regression**: `pytest tests/regression/test_laue_tof.py`.
