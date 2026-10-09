# xia2N

Fork of xia2 on branch `neutron`, adding `xia2.laue_tof` for neutron Laue
time-of-flight data. Design doc: `xia2_laue_tof_plan.md` (its Status section is
the source of truth for what is done and what is next).

## Writing style: be brief

Applies to comments, docstrings, phil `.help`, the plan, and commit messages.

- **One or two lines per comment.** A comment earns its place by saying
  something the code cannot: a constraint, a non-obvious why, an upstream bug.
  Never restate what the next line does.
- **Docstrings: one line.** Add a second short paragraph only for a real
  gotcha. Skip them entirely on short private helpers whose name says it.
- **No narrative.** Not how it was found, not what was tried first, not which
  session changed it. Git history holds that.
- **Measurements go in the plan, not the code.** A comment may carry one number
  if it is the whole justification; the rest belongs in the plan, once.
- **Say it once.** If a constraint is explained in the plan or in a phil
  `.help`, the code points at it rather than repeating it.
- Prefer deleting a comment to shortening it. Code and names first.

## Environment

`.venv/` is the interpreter; DIALS is `~/dials/dials-dev` and its dispatchers
need `export PATH=~/dials/dials-dev/conda_base/bin:$PATH` as well, because every
step shells out via `shutil.which`. See the plan's Environment section.

## Conventions

- Nothing in `Modules/Laue_TOF` may import from `Modules/SSX`; copy helpers.
- `ruff format` and `ruff check` clean; mypy clean apart from missing third-party
  stubs.
