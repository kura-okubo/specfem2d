# Porting the external-source coupling from SPECFEM2D v7.0.0 to v8.1.0

The coupling patches in this fork were written in 2019 against what was then
SPECFEM2D v7.0.0 (upstream `master`, around 2018-11-13). This branch re-applies
them to upstream **v8.1.0** (`dc4b653`, 2023-12-21, the latest release).

The fork's own history has no common ancestor with upstream — the original git
log was stripped when the fork was created — so the patch set was recovered as
`git diff ce7ea1c master` (first commit → fork tip) and re-applied hunk by hunk.

## What the coupling is

`COUPLING_IN = .true.` makes `xmeshfem2D` pick a closed band of elements
(a circle or a rectangle) and list them in `OUTPUT_FILES/externalsource.txt`;
`xspecfem2D` then reads `./extsource/EXT<iele:08d>.dat` for each of them and
prescribes that acceleration on the element's GLL nodes every time step. The
near-field solution comes from another code (HOSS); SPECFEM2D carries it to the
far field. See `EXAMPLES/Validation/note/`.

## Building and running it

### Requirements

`gfortran`, and for the examples an MPI (they use `NPROC` = 4 and 6).
`gnuplot` for the mesh check in step 2. The helper scripts under
`EXAMPLES/Validation/*/` need `numpy`, `pandas`, `matplotlib` and `scipy`;
they were last run against pandas 3.0. The test suite needs none of that.

### Build

**Linux:**

    ./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi
    make -j6 all

**macOS with Homebrew:** `configure` does not look inside the Homebrew prefix,
so the line above stops at

    checking for mpi.h... no
    configure: error: MPI header not found; try setting MPI_INC.

Pass `MPI_INC` — this is the form to use on a Mac:

    ./configure FC=gfortran CC=gcc MPIFC=mpif90 --with-mpi MPI_INC=$(brew --prefix open-mpi)/include
    make -j6 all

and export it for `make tests` as well, because the test scripts run
`configure` themselves:

    MPI_INC=$(brew --prefix open-mpi)/include make tests

(`$(brew --prefix open-mpi)/include` is `/opt/homebrew/include` on Apple
silicon, `/usr/local/include` on Intel.)

Leaving out `--with-mpi` gives a serial build, which also works but cannot run
the examples as they are shipped.

If a configure run fails, `config.log` records what it was actually given —
`grep MPI_INC config.log` shows whether the variable arrived.

### Check it works

    make tests          # on a Mac: MPI_INC=$(brew --prefix open-mpi)/include make tests

Four scripts, about 75 s: a serial build, an MPI build, an out-of-source MPI
build, and the coupling run end to end. All four have to report `[  OK  ]`.

### Run the FullSpace validation

The workflow of `EXAMPLES/Validation/note/Note_of_Specfem2DCoupling.ipynb`, at
full resolution (its six steps, with the mesh check split out as its own). Two runs of about 4 minutes each on 6 ranks, and roughly 570 MB
of intermediate files (all git-ignored).

    cd EXAMPLES/Validation/FullSpace
    rm -rf OUTPUT_FILES OUTPUT_FILES_grid OUTPUT_FILES_P OUTPUT_FILES_C extsource

`clean.sh` in the example does *not* remove those; use the line above.

**1. Mesher — pick the coupling elements.** `COUPLING_IN` is already `.true.`
in the shipped `DATA/Par_file`.

    sh run_xmeshfem2D.sh

**2. Check the band is closed.** This is the step to get right: a gap leaks.

    open OUTPUT_FILES_grid/gridfile.ps          # xdg-open on Linux
    grep -vc '^ *#' OUTPUT_FILES_grid/externalsource.txt    # -> 272

The red ring drawn on the mesh must be unbroken. If it is not, adjust `dR_ext`
(the element size usually works).

**3. Receivers on the coupling elements.**

    cd make_extsource && python make_stationfile_at_extsource.py && cd ..

**4. Point-source run — the reference.**

    sed -i '' -E 's/^(COUPLING_IN[[:space:]]*=[[:space:]]*)\.true\./\1.false./' DATA/Par_file
    sed -i '' -E 's/^(factor[[:space:]]*=[[:space:]]*)[^#]*/\1 1.0d10        /' DATA/SOURCE
    sh run_validation_pointsource.sh

**5. Turn those seismograms into the injected source.**

    cd make_extsource && python make_externalsourcefile.py && cd ..

**6. Coupling run, built-in source switched off.**

    sed -i '' -E 's/^(COUPLING_IN[[:space:]]*=[[:space:]]*)\.false\./\1.true./' DATA/Par_file
    sed -i '' -E 's/^(factor[[:space:]]*=[[:space:]]*)[^#]*/\1 0.0           /' DATA/SOURCE
    sh run_validation_extcouple.sh

**7. Compare.**

    cd plot_result && python plot_waveform_comparison.py

Black (point source) and red (coupling) should lie on top of each other at every
azimuth. `HalfSpace` is the same sequence in its own directory; it is longer
(`NSTEP` = 7000 on 4 ranks) and leaves about 1.4 GB.

The `sed -i ''` above is the BSD form, for macOS. On Linux drop the `''` and
write `sed -i -E ...`.

## Ported as-is

| file | what it does |
|---|---|
| `src/meshfem2D/determine_external_source_elements.f90` | new: picks the coupling elements, writes `externalsource.txt` |
| `src/specfem2D/compute_ext_source.F90` | new: reads the `EXT*.dat` traces and injects them |
| `src/meshfem2D/meshfem2D.F90` | calls the picker |
| `src/meshfem2D/part_unstruct.F90` | dumps `glob2loc_tableNNNNN.bin` (global element id → rank, local id) |
| `src/meshfem2D/save_databases.f90` | writes `COUPLING_IN` into the database |
| `src/meshfem2D/save_gnuplot_file.f90` | overlays the coupling elements on `gridfile.ps` |
| `src/shared/read_parameter_file.F90` | reads the `COUPLING_*` / `rec_*` parameters |
| `src/shared/shared_par.F90` | declares them, plus `iele` / `extsource` |
| `src/specfem2D/compute_forces_viscoelastic_calling_routine.F90` | calls `add_ext_source` after the mass-matrix multiply |
| `src/specfem2D/prepare_timerun.F90` | calls `read_ext_source_num` |
| `src/specfem2D/read_mesh_databases.F90` | reads `COUPLING_IN` back |
| `src/specfem2D/specfem2D_par.f90` | `glob2loc_table`, the MPI-interface tables |
| `setup/constants.h.in` | `EXT_SOURCE_NUM_MAX`, `EXT_SOURCE_TRACE_MAX` |
| `src/*/rules.mk` | build the two new objects |

## What had to change, and why

1. **`read_parameter_file.F90` moved** from `src/meshfem2D/` to `src/shared/`
   (v8 shares it between mesher and solver). The parameter block went with it.

2. **`nx` / `nz` → `nx_elem_internal` / `nz_elem_internal`, `ngnod` → `NGNOD`**
   in `determine_external_source_elements.f90`. Same quantities, renamed in v8.

3. **The MPI build macro is now `WITH_MPI`, not `USE_MPI`.** The `#ifdef USE_MPI`
   around `call smooth_MPI_interface` would have silently stopped compiling in —
   the injection would have been left un-averaged across MPI interfaces with no
   error message. This is the one change that would have been a real bug rather
   than a build failure.

4. **`module my_mpi_communicator` → `module my_mpi`**, which already re-exports
   `use mpi`, so the paired `use my_mpi_communicator` / `use mpi` collapsed to one.

5. **The MPI-only half of `compute_ext_source.F90` is now guarded** by
   `#ifdef WITH_MPI`, and so is the call to `setup_iglob_interface`. In v7 those
   routines had a bare `use mpi` that broke a serial (non-MPI) build. Both builds
   compile now.

6. **`read_ext_source_num` distributes with `bcast_all_*` instead of
   `send_i` / `recv_i`.** v8 only provides the point-to-point wrappers in an MPI
   build, so the serial link failed. A broadcast is equivalent and is O(1) calls
   instead of O(NPROC × number_of_extsource).

7. **The coupling parameters are optional.** v7 made all eleven mandatory in
   every `Par_file`, which would break every stock example. They now default
   (`COUPLING_IN = .false.`, the rest 0) and only the ones behind
   `if (COUPLING_IN)` are required. Upstream examples run untouched.

8. **`read_mesh_databases.F90` validates rather than overwrites.** v8 reads each
   database flag into a local and warns if it disagrees with the `Par_file`;
   `COUPLING_IN` now follows that pattern.

9. **`construct_glob2loc_elmnts`**: dropped the `write(*,*) num_loc(1)` debug
   lines, which index out of bounds when `NPROC = 1`.

## Not ported

Three v7 hunks were personal debugging aids, not part of the coupling. They are
kept as patches under `optional_patches/` rather than carried into the branch:

* `src/specfem2D/check_grid.F90` — forced the PostScript grid plot on. v8 has
  `output_postscript_snapshot` in the `Par_file`; set that instead.
  (`optional_patches/force_grid_postscript.patch`)
* `src/specfem2D/write_wavefield_dumps.F90` — dumped element corners only
  instead of all GLL points. (`optional_patches/dump_corner_points_only.patch`)
* `src/specfem2D/define_external_model_from_tomo_file.f90` — `write(IMAIN,...)`
  per GLL point; dropped outright.

`src/shared/exit_mpi.F90`'s `ERROR stop 1` is now unnecessary: v8 already ends
`exit_MPI` with `stop 30`.

## Example updates

The `EXAMPLES/Validation/*` inputs were written for v7 and needed the v8 format:

* `DATA/Par_file`: `ngnod` → `NGNOD`, `partitioning_method` → `PARTITIONING_TYPE`,
  `NSTEP_BETWEEN_*` → `NTSTEP_BETWEEN_*`, `subsamp_seismos` →
  `NTSTEP_BETWEEN_OUTPUT_SAMPLE`; added `noise_source_time_function_type`,
  `write_moving_sources_database`, `APPROXIMATE_HESS_KL`.
* `DATA/SOURCE`: v8 wants 16 lines per source, v7 had 14 — added `vx`, `vz`
  (moving sources). Without them the mesher aborts with
  "invalid number of (non-blank and non-comment) lines per source".
* The helper scripts: `delim_whitespace=` was removed in pandas 3.0, so they now
  use `sep=r'\s+'`; `make_stationfile_at_extsource.py` in `FullSpace_reccoupling`
  had `temp_df.falues`.
* `run_xmeshfem2D.sh` now runs `gnuplot OUTPUT_FILES/plot_gridfile.gnu` from the
  example root. v8.1.0's `save_gnuplot_file.f90` writes `./OUTPUT_FILES/...` into
  that script, so the old `cd ./OUTPUT_FILES && gnuplot plot_gridfile.gnu` made
  every path in it wrong. (v7 had patched the paths in the Fortran instead;
  fixing the caller keeps the source diff smaller.)

## Validation

Both examples from the note were run end to end (gfortran 15.2, Open MPI 5.0.9,
macOS arm64), following the six steps in `note/Note_of_Specfem2DCoupling.ipynb`.

| | coupling elements | note's v7 value | result |
|---|---|---|---|
| `FullSpace` (circle, R = 5 km, NPROC = 6) | 272 | — | coupling reproduces the point source |
| `HalfSpace` (circle, R = 20 km, PML, NPROC = 4) | 510 | 510 | ditto |

`HalfSpace` finding exactly the 510 coupling elements the note records is a
direct check that the mesher-side port is faithful.

Agreement between the point-source reference and the coupling run, over the 12
validation receivers at R = 20 km, after the 10 Hz low-pass the note's figure
uses: **median 9.6 %, worst 10.1 % L2 misfit** for `FullSpace`. Unfiltered the
misfit is 12-53 %, which is the high-frequency numerical oscillation the note
already describes ("we may need grid convergence analysis as there is a small
error in amplitude and phase associated with coupling model").

`gridfile.ps` comes out with the coupling elements drawn in red on top of the
mesh, closed, exactly as `note/fig/gridfile.png` shows.

## Azimuth axis

The validation traces are now labelled by **azimuth measured clockwise from
north**, with north = `+z` and east = `+x`, matching the map view the coupled
application uses (HOSS `x`-`y` → SPECFEM2D `x`-`z`). 0° = N, 90° = E,
180° = S, 270° = W.

The receiver positions did not move — `theta` is still only the parameter that
generates the ring, and `DATA/STATIONS` comes out byte for byte identical, so
runs already on disk stay valid. What this fixed along the way: the
FullSpace-type `plot_waveform_comparison.py` built its ring at `theta` while
`make_stationfile_at_extsource.py` wrote the receivers at `theta - 90°`, so the
old axis was 90° away from the receiver it labelled. (`HalfSpace` already had
the two scripts in agreement.)

Note for `HalfSpace`: that model is a vertical section with a free surface at
`z = 0`, so "azimuth from north" there is really the angle from `+z`, i.e. from
straight up. The label is consistent with `FullSpace`; just do not read it as a
map azimuth.

## EXAMPLES

Everything under `EXAMPLES/` except `Validation/` was deleted (48 entries,
~340 MB tracked). Consequences worth knowing:

* The CI configurations still name examples that are gone:
  `.github/workflows/CI.yml` (8 paths), `.github/scripts/run_tests.sh` (4),
  `.travis/run_tests.sh` (21), `.azure-pipelines.yml` (2), `.travis.yml` (1),
  and both scripts under `tests/examples/`. They were left untouched; they will
  fail if anything ever runs them.
* `EXAMPLES/process_DATA_Par_files_to_update_their_format_when_new_parameters_are_added.py`
  was kept: it is upstream's own tool for exactly the Par_file upgrade done here.

## CI and tests

`make tests` is now the whole test suite, and it is green: four scripts,
about 75 s end to end.

    tests/compilations/0.configure.serial_make.sh     serial build
    tests/compilations/1.configure.parallel_make.sh   MPI build
    tests/examples/0.configure.parallel_make.sh       MPI build, out of source
    tests/examples/1.run_coupling_validation.sh       the coupling, end to end

`tests/examples/1.run_coupling_validation.sh` replaces
`1.run_simple_topography_and_also_a_simple_fluid_layer.sh`, whose example is
gone. It runs the six-step workflow of the note on a shrunk version of
`EXAMPLES/Validation/FullSpace` — 50 km box at 500 m elements, 14 s of
propagation, 2 MPI ranks, 128 coupling elements — and requires the coupling run
to correlate with the point-source run at receivers outside the coupling
surface. Roughly 13 s.

There is deliberately **no stored reference solution**: the test compares the
coupling against a point-source run it computes in the same job, which is the
statement the method has to satisfy. Nothing to regenerate when the numerics
shift slightly, and nothing large to track.

Measured correlations are 0.983-0.991; the threshold is 0.8, because the
injection leaves a high-frequency oscillation that the note's own figures
low-pass away before comparing. The test was checked against a sabotaged run
(EXT traces zeroed): correlation drops to 0.0 and it exits 1, so it is not a
rubber stamp.

Two dependencies were deliberately avoided so the test needs nothing beyond a
POSIX shell: the two `pandas` helper scripts are replaced by `awk`, and `sed -i`
is wrapped in a `sedi()` that works with both GNU and BSD `sed` (upstream's
tests assume GNU `sed` and rely on CI installing `gnu-sed` on macOS).

`.github/workflows/CI.yml` keeps `changesCheck`, `macosCheck`, `linuxCheck` and
the `make tests` job. The eleven `linuxTest_1..11` jobs were removed: each ran
one of the deleted examples through `.github/scripts/run_tests.sh`, which needs
a `run_this_example.sh` and a `REF_SEIS/` that no Validation example has. That
script is left in place, unused, in case examples ever come back.

Four more things had to change before the workflow would go green, none of them
caused by this port:

* `changesCheck` fetched `github.event.before` unconditionally. On the first
  push of a new branch that is the null SHA; the fetch failed, and since every
  other job needs `changesCheck`, the whole workflow was skipped. It now falls
  back to the files of the pushed commit.
* `macosCheck` runs `make tests`, which runs `configure --with-mpi` itself, and
  `configure` does not search the Homebrew prefix — the same `MPI header not
  found` seen when building by hand. `MPI_INC` is now exported in the install
  step. This only started failing when `macos-latest` moved to Apple silicon:
  on the old Intel runners the prefix was `/usr/local`, which is searched by
  default.
* The `ubuntu-20.04` runner image has been retired, so jobs asking for it queue
  forever instead of failing. The `linuxCheck` matrix is now
  `[ubuntu-latest, ubuntu-22.04]`, and `linuxCheck-Intel`, which pinned
  `ubuntu-20.04` and installed oneAPI 2023.2.2, was removed — a long, fragile
  step for a compiler this fork does not use.
  `git checkout v8.1.0 -- .github/workflows/CI.yml` has it back.
* `actions/checkout` moved from v3 (Node 16, retired) to v4.

`.travis.yml`, `.travis/`, `.azure-pipelines.yml` and `.azure-pipelines/` were
deleted. They were not connected to this fork, and between them they named 23
examples that no longer exist. `git checkout v8.1.0 -- <path>` restores any of
them.

`README.md` regained the fork's own header, now also pointing at this file, and
its CI badge points at `kura-okubo/specfem2d` instead of upstream, so it reports
this fork's status. The dead Travis and Azure badges were dropped.
