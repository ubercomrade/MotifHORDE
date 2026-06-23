# SiteGA pybind11 integration plan

## Goal

Make SiteGA a Python extension module built with the `motifhorde` package and
use it directly from `SitegaDiscoveryTool`. SiteGA must no longer depend on the
external `andy05cell.exe` executable or a `PATH` lookup.

## Current state

- The main project uses `setuptools.build_meta` with a `src` layout.
- SiteGA sources are in `src/sitega/`:
  - `andy05cell.cpp`
  - `bindings.cpp`
  - `sitega_train.h`
  - `setup.py`
- `bindings.cpp` exposes a `sitega` Python module with:
  - `sitega.TrainParams`
  - `sitega.train(...) -> (rc, mat_path, loc_path, best_fit)`
- `SitegaDiscoveryTool` still calls `run_sitega(...)`, which shells out to
  `andy05cell.exe`.
- The CLI still checks for `andy05cell.exe`; this check should be removed.

## Implementation plan

### 1. Move SiteGA extension build into the main package

Keep one packaging entry point for the repository.

Changes:

- Add `pybind11` to `[build-system].requires` in `pyproject.toml`.
- Add a root-level `setup.py` or equivalent setuptools hook only for extension
  module configuration.
- Build extension module name: `sitega`.
- Build sources:
  - `src/sitega/bindings.cpp`
  - `src/sitega/andy05cell.cpp`
- Include directory:
  - `src/sitega`
- Use C++17.

Rationale:

- `pyproject.toml` remains the source of project metadata.
- The root build config makes `sitega` part of the normal `uv sync`,
  `pip install .`, and wheel build flow.
- `src/sitega/setup.py` becomes redundant and should be removed after the root
  build is working.

### 2. Handle compiler and platform differences explicitly

Do not hard-code Linux/GCC-only flags without a fallback.

Preferred build behavior:

- Always compile with C++17.
- Use conservative optimization first, for example `-O3` on GCC/Clang and
  `/O2` on MSVC.
- Avoid enabling `-ffast-math` by default unless there is a measured need. It
  can change floating-point behavior and conflicts with the project preference
  for correctness over clever optimization.
- Enable OpenMP only when the platform/compiler can support it.
- Compile without OpenMP as a supported fallback. The C++ code already guards
  OpenMP-specific calls with `_OPENMP`.

Platform notes:

- Linux with GCC:
  - likely flags: `-O3`, `-fopenmp`
  - link flags: `-fopenmp`
  - this should be the primary tested path in the current environment.
- Linux with Clang:
  - `-fopenmp` may require `libomp` and extra include/library paths.
  - if detection is not implemented, fall back to no OpenMP instead of failing
    the build.
- macOS with Apple Clang:
  - OpenMP is usually unavailable by default.
  - require `libomp` only if OpenMP support is explicitly enabled.
  - default fallback should be a serial build.
- Windows with MSVC:
  - use `/std:c++17`, `/O2`, optionally `/openmp`.
  - avoid POSIX-only assumptions in build configuration.
  - runtime behavior still needs validation because the C++ code uses fixed-size
    char buffers and path separators are normalized only partly in bindings.

Pragmatic first implementation:

- Add a small build helper that chooses flags by `sys.platform` and compiler
  type.
- Default to OpenMP on Linux/GCC.
- Fall back to serial build on unsupported compilers instead of blocking
  installation.
- Keep this logic local to packaging; do not introduce a larger build framework
  unless setuptools becomes insufficient.

### 3. Replace subprocess SiteGA execution

Update `src/motifhorde/discovery.py`.

Changes:

- Remove `from .execute import run_sitega`.
- Import `sitega` lazily inside the SiteGA execution path, not at module import
  time.
- Call `sitega.train(...)` directly with full FASTA paths.
- Pass:
  - `fg_path=foreground`
  - `bg_path=background`
  - `max_lpd=6`
  - `motif_len=length`
  - `size=number_of_lpd`
  - `olig_bg=6`
  - `infc=1`
  - `out_path=output_dir + os.sep`
  - `max_peak_len=5000`
  - `log_file="sitega.log"`
- Check the returned `rc`.
- Raise `RuntimeError` with a clear message if `rc != 0`.
- Prefer reading the returned `mat_path`.
- Keep a fallback glob for `train.fa_mat*` only if the binding returns an empty
  path despite success.

Rationale:

- The discovery wrapper should not copy FASTA files just to satisfy the old CLI
  calling convention.
- The pybind11 binding already splits full paths into directory and filename.
- Reading the returned `mat_path` makes output handling explicit and less
  dependent on filename conventions.

### 4. Remove old executable runner

Update `src/motifhorde/execute.py`.

Changes:

- Remove `run_sitega(...)`.
- Remove imports that become unused after deleting it.

Rationale:

- Keeping both subprocess and pybind11 paths creates hidden behavior drift.
- SiteGA should have one execution path.

### 5. Remove CLI dependency check for `andy05cell.exe`

Update `src/motifhorde/cli.py`.

Changes:

- Delete the `sitega` branch that checks `shutil.which("andy05cell.exe")`.
- Do not replace it with a CLI-level binary check.

Rationale:

- SiteGA is now a Python extension module.
- If installation is broken, discovery should fail with a Python import/build
  error, not a misleading missing executable error.
- This keeps dependency validation consistent with packaged Python code.

### 6. Update docs

Update `README.md`.

Changes:

- Replace the note that `sitega` requires `andy05cell.exe` in `PATH`.
- Document that SiteGA is built as a pybind11 extension during installation.
- Mention build requirements:
  - C++17 compiler
  - `pybind11` as a build dependency
  - optional OpenMP support
- Document that non-OpenMP builds are supported but may be slower.

### 7. Update tests

Add or update tests around the changed contract.

Unit tests:

- Add a `SitegaDiscoveryTool` test that monkeypatches a fake `sitega` module.
- Fake `sitega.train(...)` should return `(0, mat_path, loc_path, best_fit)`.
- Create a minimal `.mat` fixture compatible with `mimosa.read_model(...,
  "sitega")`.
- Assert that returned motifs are named `Sitega-1`, filtered by expected length,
  and read from the returned `mat_path`.

CLI tests:

- Add a test that `check_external_dependencies` or the equivalent CLI path does
  not require `andy05cell.exe` for `sitega`.
- Keep existing `--lpd` parameter parsing tests.

Build smoke tests:

- Add a lightweight optional smoke check:
  - `uv run python -c "import sitega; print(sitega.__doc__)"`
- Keep it separate from normal unit tests if extension compilation is expensive
  or platform-sensitive.

### 8. Update lock/environment files if needed

Changes:

- Regenerate `uv.lock` after adding `pybind11` to build requirements if uv
  changes the lock file.
- Add `conda-forge::openmp` to `environment.yml` so conda-based installs have
  the OpenMP runtime available for the SiteGA extension.
- `environment.yml` does not need a runtime `pybind11` dependency if builds are
  performed through PEP 517 isolation.
- Keep compiler/OpenMP notes in docs because non-conda installs still depend on
  the platform toolchain and runtime libraries.

## Verification plan

Run local checks in this order:

```bash
uv sync
uv run python -c "import sitega; print(sitega.__doc__)"
uv run pytest tests/test_cli_tools.py
uv run pytest tests/test_sitega_discovery.py
uv run ruff check .
```

Optional full pipeline check:

```bash
HORDEMOTIFS_RUN_FULLRUN=1 uv run pytest tests/test_external_fullrun.py -k sitega
```

## Main risks and mitigations

- OpenMP flags are compiler-specific.
  Mitigation: add `conda-forge::openmp` for conda environments, detect
  platform/compiler in the build, and fall back to serial build when OpenMP is
  unavailable.

- macOS may not have OpenMP available.
  Mitigation: do not require OpenMP by default on macOS; document optional
  `libomp` setup later if needed.

- Windows build may expose path or compiler assumptions in the adapted C++ code.
  Mitigation: keep the first implementation portable at the setuptools flag
  level, then validate Windows separately before claiming support.

- Fixed-size C buffers in `andy05cell.cpp` may fail on very long paths.
  Mitigation: use short output paths in tests and document this as a known C++
  limitation if not fixed now.

- `sitega` top-level module name could conflict with an unrelated installed
  package named `sitega`.
  Mitigation: acceptable for now because the binding is intentionally imported
  as `sitega`; revisit only if an actual conflict appears.

- Build isolation may hide missing compiler errors until install time.
  Mitigation: keep the build dependency list explicit and add the import smoke
  check to verification.

## Definition of done

- `sitega` imports after installing the project.
- `SitegaDiscoveryTool` calls `sitega.train(...)` directly.
- No code path calls `andy05cell.exe`.
- CLI no longer checks for `andy05cell.exe`.
- `run_sitega(...)` is removed.
- Documentation reflects pybind11 extension build requirements.
- Unit tests cover the new Python-module execution path.
- Ruff and targeted pytest checks pass.
