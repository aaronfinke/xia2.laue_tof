from __future__ import annotations

phil_str = """
image = None
  .type = str
  .multiple = True
  .help = "Path to an image file. Give one per exposure of the run; they are"
          "stationary measurements of one crystal at different settings, and"
          "are processed together."
  .expert_level = 0
template = None
  .type = str
  .multiple = True
  .help = "The image sequence template"
  .expert_level = 0
directory = None
  .type = str
  .multiple = True
  .help = "A directory with images"
  .expert_level = 0
vanadium_run = None
  .type = path
  .help = "Incident run (e.g. vanadium), used to normalise the intensities by the"
          "incident spectrum during integration. Optional, but recommended."
  .expert_level = 0
empty_run = None
  .type = path
  .help = "Empty instrument run, used to correct the incident run for background."
          "Only used if vanadium_run is also given."
  .expert_level = 0
mask = None
  .type = path
  .help = "A mask to use for spotfinding and integration"
  .expert_level = 0
detector_calibration = None
  .type = path
  .help = "An ISAW/Mantid .DetCal file giving refined per-panel positions."
          "Applied straight after import, so indexing, refinement and"
          "integration all use the calibrated geometry. Without it the nominal"
          "geometry can be wrong by several pixels per panel, which refinement"
          "absorbs into the crystal orientation rather than correcting, leaving"
          "per-panel residuals that no mask or profile model can recover from."
  .expert_level = 0
space_group = None
  .type = space_group
  .help = "Space group to be used for indexing and integration."
  .expert_level = 0
d_min = None
  .type = float
  .help = "High resolution limit of the data that is kept, in Angstrom."
          "Applied to the integrated reflections by their own d spacing, so it"
          "reaches the SHELX export, the MTZ and the scaling without"
          "re-integrating. Not a spotfinding cut: a DIALS one is a pixel mask,"
          "and a Laue-TOF pixel has no single resolution. Where to cut is a"
          "merging decision - prefer to leave this unset and take it from"
          "CC-half after scaling."
  .expert_level = 0
d_max = None
  .type = float
  .help = "Low resolution limit of the data that is kept, in Angstrom. As"
          "d_min, applied to the integrated reflections rather than to the"
          "spot list."
  .expert_level = 0
multiprocessing {
  nproc = Auto
    .type = int
    .expert_level = 1
}
dials_import.phil = None
  .type = path
  .help = "Phil file to use for dials.import. Parameters defined in the"
          "xia2.laue_tof phil scope will take precedent over identical options"
          "defined in the phil file."
  .expert_level = 2
"""

spotfinding_phil_str = """
spotfinding {
  min_spot_size = 15
    .type = int
    .help = "The minimum spot size to allow in spotfinding."
    .expert_level = 0
  max_spot_size = None
    .type = int
    .help = "The maximum spot size to allow in spotfinding."
    .expert_level = 1
  threshold_algorithm = dispersion *dispersion_extended radial_profile
    .type = choice
    .expert_level = 1
  max_strong = 35000
    .type = int
    .help = "If more than this many strong spots are found, stop rather than"
            "proceed to indexing."
    .expert_level = 2
  phil = None
    .type = path
    .help = "Phil options file to use for spotfinding. Parameters defined in"
            "the xia2.laue_tof phil scope will take precedent over identical"
            "options defined in the phil file."
    .expert_level = 2
}
"""

indexing_phil_str = """
indexing {
  unit_cell = None
    .type = unit_cell
    .expert_level = 0
  method = None
    .type = str
    .help = "Indexing method to use. If None, each method given in"
            "indexing.ladder is tried in turn until one succeeds."
    .expert_level = 0
  ladder = fft3d real_space_grid_search fft1d
    .type = strings
    .help = "Indexing methods to try, in order, when indexing.method is not set."
            "fft1d is the method that succeeds on sparse, narrow-wedge data from"
            "distant detector layouts, where fft3d cannot populate its grid."
    .expert_level = 1
  max_cell = None
    .type = float
    .help = "Fixed maximum cell length (Angstrom) for the basis vector search."
            "If None, it is set to max_cell_multiplier times the longest edge of"
            "indexing.unit_cell when a cell is given, and to max_cell_fallback"
            "when one is not. The DIALS nearest-neighbour estimate is never"
            "used: it is unreliable on sparse Laue-TOF spot lists, and when it"
            "comes out too small it does not fail, it returns a truncated cell"
            "that passes every acceptance test."
    .expert_level = 1
  max_cell_multiplier = 1.3
    .type = float(value_min=1.0)
    .expert_level = 2
  max_cell_fallback = 150
    .type = float
    .help = "The max_cell used when neither indexing.max_cell nor"
            "indexing.unit_cell is given. A fixed value beats the"
            "nearest-neighbour estimate, which is unreliable on sparse"
            "Laue-TOF spot lists: on NMX rubredoxin it returned 23.3 A for a"
            "44.1 A axis and halved the cell, while still passing every"
            "acceptance check. Raise it for a longer edge than this, or better,"
            "give indexing.unit_cell and make it moot."
    .expert_level = 1
  d_min_start = 4.0
    .type = float
    .help = "Resolution limit for the first macrocycle of indexing refinement."
    .expert_level = 1
  min_spots = 10
    .type = int
    .help = "Only attempt to index an exposure with at least this many strong spots."
    .expert_level = 1
  min_indexed = 50
    .type = int
    .help = "Reject an indexing solution with fewer than this many indexed"
            "reflections, and move on to the next method in the ladder. This is"
            "a floor: min_indexed_fraction of the strong spots is used when that"
            "is the larger number."
    .expert_level = 1
  min_indexed_fraction = 0.25
    .type = float(value_min=0, value_max=1)
    .help = "Reject an indexing solution that indexes a smaller fraction of the"
            "strong spots than this. A fraction travels between datasets in a"
            "way an absolute count does not: a good MANDI solution indexes 65%"
            "of 227 spots, which an absolute cut of 250 would have thrown away,"
            "while a bad solution on 2300 NMX spots can index 300 of them."
    .expert_level = 1
  max_rmsd_px = 5.0
    .type = float
    .help = "Reject an indexing solution whose positional residual is worse"
            "than this, in pixels. Indexing plenty of reflections on a cell"
            "of the right size is not enough: fft3d has taken 39% of a MANDI"
            "spot list at 6.1 px where another method fitted it at 3.2 px."
            "The median, not the mean, which is dominated by the tail a good"
            "solution can still leave. None turns the test off."
    .expert_level = 1
  target_rmsd_px = 2.0
    .type = float
    .help = "An accepted solution this good stops the ladder. Otherwise every"
            "method is tried and the best by RMSD kept, the first acceptable"
            "answer not necessarily being the best fit. None makes the first"
            "one win."
    .expert_level = 1
  seed_from_first = True
    .type = bool
    .help = "After every method in the ladder has been tried, index whatever is"
            "still unsolved from the crystal model of the best solution found."
            "All the exposures are measurements of the same crystal, so this"
            "both enforces a consistent setting and recovers exposures that"
            "will not index on their own."
    .expert_level = 0
  outlier {
    algorithm = auto *mcd tukey sauter_poon null
      .type = choice
      .help = "Outlier rejection algorithm used in indexing refinement. The"
              "residuals of sparse Laue-TOF data carry a heavy tail that a"
              "tight cut is needed to remove."
      .expert_level = 1
    iqr_multiplier = 0.8
      .type = float
      .expert_level = 2
    separate_panels = True
      .type = bool
      .expert_level = 2
  }
  phil = None
    .type = path
    .help = "Phil options file to use for indexing. Parameters defined in the"
            "xia2.laue_tof phil scope will take precedent over identical options"
            "defined in the phil file."
    .expert_level = 2
}
"""

integration_phil_str = """
integration {
  method = profile_1d_ibix
    .type = str
    .help = "Integration method, as named by dials.tof_integrate, which is"
            "asked at startup which ones it offers rather than the list being"
            "repeated here. Summation intensities are written whichever method"
            "is chosen; a profile method adds intensity.prf alongside them."
            "Profile fitting is not robust on all data - it can abort the whole"
            "run on a bad shoebox - so see fallback_method."
    .expert_level = 0
  fallback_method = summation
    .type = str
    .help = "What to run when a profile method fails outright, or none to let"
            "the failure stop the run. The failure and the fallback are"
            "logged. The whole run is integrated in one call, so the fallback"
            "applies to every exposure at once."
    .expert_level = 1
  min_profile_fraction = 0.5
    .type = float(value_min=0, value_max=1)
    .help = "The fraction of reflections a profile method has to fit for its"
            "intensities to be usable. A profile run can finish successfully"
            "having fitted almost nothing - profile_3d_gutmann fitted 13 of"
            "1359 reflections on NMX test data - which is worse than useless,"
            "because the profile column then exists but covers a handful of"
            "reflections. Below this fraction the summation intensities are"
            "used instead."
    .expert_level = 1
  integration_type = *observed calculated
    .type = choice
    .help = "observed integrates only the spots found by spotfinding."
            "calculated integrates everything the geometry predicts out to"
            "calculated_d_min - two to three orders of magnitude more, most of"
            "it weak, which is the point: where to cut is a merging decision."
            "A calculated run makes the observed pass first, to choose the"
            "method and mask cheaply, and keeps it as a checkpoint."
    .expert_level = 1
  calculated_d_min = None
    .type = float
    .help = "Resolution limit used when integration_type=calculated. If None,"
            "the limit the detector and the wavelength band can reach is used"
            "- d = lambda / (2 sin theta) at the shortest wavelength and the"
            "largest scattering angle any pixel sees - since nothing beyond it"
            "can have been recorded."
    .expert_level = 1
  max_predicted = 250000
    .type = int(value_min=1)
    .help = "Refuse a calculated pass that would integrate more reflections"
            "than this. Prediction takes seconds and integration takes hours,"
            "so the count is checked first: a large cell with an optimistic"
            "d_min can predict millions. The message says which d_min would"
            "bring the count under the guard."
    .expert_level = 1
  background_model = constant2d constant3d linear2d *linear3d
    .type = choice
    .expert_level = 1
  mask = ellipse seed_skewness *both
    .type = choice
    .help = "Foreground/background mask method. both integrates twice and"
            "keeps the better by mask_metric. The masks partition the shoebox"
            "differently, so which is better is data dependent: seed_skewness"
            "takes in more of the peak tails and gives higher I/sigma, at"
            "several times the run time. One call covers the run, so every"
            "batch is integrated alike."
    .expert_level = 1
  mask_metric = *i_over_sigma n_integrated
    .type = choice
    .help = "How to choose between the two masks when mask=both. Neither is a"
            "substitute for merging statistics, which are only available after"
            "scaling: treat the choice as a starting point and check it there."
    .expert_level = 2
  ellipse_mask_scale = 1.0
    .type = float(value_min=0.5)
    .expert_level = 2
  wavelength_range = None
    .type = floats(size=2)
    .help = "Reflections outside this wavelength range (Angstrom) are not integrated."
    .expert_level = 1
  bbox_tof_padding = 2
    .type = int
    .expert_level = 2
  bbox_xy_padding = 1
    .type = int
    .expert_level = 2
  absorption {
    enabled = False
      .type = bool
      .help = "Apply the spherical absorption correction. Requires the target"
              "spectrum parameters below to be set."
      .expert_level = 1
    sample_number_density = None
      .type = float
      .help = "Sample number density (atoms/A^3) of the crystal."
      .expert_level = 1
    sample_radius = None
      .type = float
      .help = "Sample radius (mm) of the crystal."
      .expert_level = 1
    scattering_x_section = None
      .type = float
      .help = "Sample scattering cross section (barns) of the crystal."
      .expert_level = 1
    absorption_x_section = None
      .type = float
      .help = "Sample absorption cross section (barns) of the crystal."
      .expert_level = 1
  }
  phil = None
    .type = path
    .help = "Phil options file to use for integration. Parameters defined in the"
            "xia2.laue_tof phil scope will take precedent over identical options"
            "defined in the phil file."
    .expert_level = 2
}
"""

output_phil_str = """
output {
  intensity = *auto sum profile
    .type = choice
    .help = "Which intensities to export. auto uses the profile-fitted"
            "intensities when they are present and cover at least"
            "integration.min_profile_fraction of the reflections, and the"
            "summation intensities otherwise. dials.export refuses its own auto"
            "when both columns exist, so one is always chosen here explicitly."
            "The unmerged MTZ carries summation intensities as well when"
            "profile is chosen."
    .expert_level = 0
  composition = CH
    .type = str
    .help = "Chemical composition of the asymmetric unit, written into the SHELX"
            ".ins file that accompanies the exported intensities."
    .expert_level = 1
  phil = None
    .type = path
    .help = "Phil options file to use for the export. Parameters defined in the"
            "xia2.laue_tof phil scope will take precedent over identical options"
            "defined in the phil file."
    .expert_level = 2
}
"""

scaling_phil_str = """
scaling {
  lawless {
    executable = None
      .type = path
      .help = "Path to lawless (the Laue branch of aimless). Defaults to"
              "'lawless' then 'aimless' on $PATH."
      .expert_level = 0
    probe = *neutron xray
      .type = choice
      .expert_level = 1
    normalisation = *gpr chebyshev none
      .type = choice
      .help = "How the incident spectrum is modelled: a Gaussian process"
              "(NORMGPR), a Chebyshev polynomial (NORMCHEBYSHEV) or not at all."
      .expert_level = 0
    lam_min = None
      .type = float
      .help = "Wavelength range the normalisation curve covers. Defaults to the"
              "range of the data."
      .expert_level = 1
    lam_max = None
      .type = float
      .expert_level = 1
    lam_ref = None
      .type = float
      .help = "Reference wavelength the curve is normalised to (NORMLAMREF)."
              "Defaults to the middle of the range."
      .expert_level = 1
    chebyshev_degree = 6
      .type = int
      .expert_level = 1
    gpr_bins = None
      .type = int
      .help = "Training bins for the Gaussian process (NORMGPRBINS). Too many"
              "gives a curve with structure finer than the physics justifies."
      .expert_level = 1
    scales = "BATCH BFACTOR ON SECONDARY 0"
      .type = str
      .help = "The scale model, as the SCALES keyword takes it. One scale and"
              "one B factor per batch, which with a handful of stationary"
              "exposures may already be over-parameterised."
      .expert_level = 0
    sdcorrection = "NOREFINE 1.0 0.0 0.0"
      .type = str
      .help = "SDCORRECTION keyword. The default leaves the sigmas alone, so"
              "that I/sigma stays comparable with the unscaled data."
      .expert_level = 1
    anomalous = False
      .type = bool
      .expert_level = 1
    resolution = None
      .type = floats(size=2)
      .help = "Resolution range to scale over. The cutoff is a merging"
              "decision, so it belongs here rather than at integration."
      .expert_level = 0
    lambda_only = False
      .type = bool
      .help = "Run the LAMBDAONLY diagnostic, which fits the wavelength curve"
              "and stops."
      .expert_level = 1
    keywords = None
      .type = str
      .multiple = True
      .help = "Extra lawless keyword lines, appended after the generated ones,"
              "so that they win."
      .expert_level = 1
  }
}
"""

workflow_phil_str = """
workflow {
  steps = *find_spots *index *refine *integrate *combine *export \
          *unmerged_mtz *lawless *report
    .type = choice(multi=True)
    .help = "Option to turn off particular steps. Multiple choices should be of"
            "the format steps=find_spots+index"
    .expert_level = 2
}
"""

full_phil_str = (
    phil_str
    + spotfinding_phil_str
    + indexing_phil_str
    + integration_phil_str
    + output_phil_str
    + scaling_phil_str
    + workflow_phil_str
)
