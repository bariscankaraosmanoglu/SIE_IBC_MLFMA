# SIE_IBC_MLFMA

A C++ electromagnetic scattering solver for targets with impedance boundary
conditions (IBC), using surface integral equations, Rao–Wilton–Glisson (RWG)
basis functions on triangular meshes, and the Multilevel Fast Multipole
Algorithm (MLFMA). Author: **Bariscan Karaosmanoglu, PhD**.

## Scope

The coupled system has **2N unknowns for N RWGs**: the first N coefficients
represent electric current scaled by the free-space impedance (`eta0 * J`),
and the second N represent magnetic current (`M`). Four stored near matrices
and two far-operator branches couple these currents. The existing IBC signs,
quadrature, impedance terms, and field normalization are retained.

The executable supports one frequency per run, plane-wave illumination with
incident-direction sweeps, GMRES or flexible GMRES, and direct scattered
far-field output. Parallelism uses shared-memory OpenMP; there is no MPI or GPU
backend. Reusable workspace caches angular tables, interpolation stencils, and
eight child-shift phases per level. Full per-RWG radiation-pattern caching is
not used.

Input names remain compatible with the related PEC solver. This executable
accepts only `Formulations 1` and `source_type 0`, the implemented IBC and
plane-wave modes. It rejects unsupported selectors explicitly. PEC CFIE weights
and PEC open-component EFIE forcing are not applied to this coupled IBC system.

## Build

Requirements:

- GNU Make and a C++ compiler supporting GNU C++17 and OpenMP.
- The compiler's OpenMP runtime (for GCC, `libgomp`).
- No external BLAS, LAPACK, or mesh library is required. Python 3 is needed only
  for the optional regression scripts.

From the repository root:

```sh
make -j2
```

This produces `./run_mlfma`; intermediate objects and generated header
dependencies are stored in `build/`. Choose the `-j` count to suit available
build memory. Build jobs are independent of the solver's runtime thread count.
Building does not run the solver or tests.

Default options preserve the existing build:

```sh
g++ -Wall -fexceptions -O3 -std=gnu++17 -fopenmp src/*.cpp -o run_mlfma
```

The Makefile compiles individual translation units and links them incrementally.
Header changes rebuild affected objects; changes to the compiler, flags, or
source list rebuild the configuration. GCC 11.4.0 and GNU Make 4.3 have been used
locally. Existing source warnings are not treated as errors.

Standard overrides are available: `CXX`, `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`,
`LDLIBS`, plus `OPENMP_FLAGS`, `BUILD_DIR`, and `TARGET`. OpenMP options are
applied during compilation and linking. Keep GNU C++17 and OpenMP enabled.
For a separate debug build:

```sh
make -j2 BUILD_DIR=build/debug TARGET=build/debug/run_mlfma \
  CXXFLAGS='-Wall -fexceptions -O0 -g -std=gnu++17'
```

`make help` lists targets. `make clean` explicitly removes the selected
configuration's objects, dependency files, configuration record, and executable;
it leaves simulation inputs, outputs, and logs intact. Supply the same
`BUILD_DIR` and `TARGET` overrides when cleaning a custom build. No cleanup runs
automatically.

### Optional compiler optimizations

The default build remains unchanged. The following optional GCC switches apply
at both compilation and linking, and changes trigger a rebuild:

| Make variable | Effect |
| --- | --- |
| `LTO=2` | Enable link-time optimization with up to two compiler jobs at the LTO stage. Default `LTO=0` disables it. |
| `TUNE=native` | Tune for the build CPU while retaining the selected instruction set. Default: no additional tuning flag. |
| `ARCH=native` | Enable the build CPU's instruction set. The executable may require that CPU's features on other systems. Default: no additional architecture flag. |
| `PGO=generate` | Build an instrumented executable for profile collection, using atomic counter updates for OpenMP workers. |
| `PGO=use` | Rebuild using collected profiles; missing or mismatched profiles are compilation errors. Default `PGO=off` disables profiling. |
| `PROFILE_DIR=path` | Profile-data location, resolved to an absolute path. Default: `$(BUILD_DIR)/profiles`. |

`ARCH` and `TUNE` also accept explicit GCC CPU names. For an initial LTO/tuning
candidate that retains the default instruction set:

```sh
make -j2 LTO=2 TUNE=native \
  BUILD_DIR=build/lto TARGET=build/lto/run_mlfma
```

For an executable intended for the build machine's CPU features:

```sh
make -j2 LTO=2 ARCH=native \
  BUILD_DIR=build/native TARGET=build/native/run_mlfma
```

The separate directories retain the default executable and build products.
LTO can increase build time and compiler memory. CPU-specific instructions,
particularly fused multiply-add, can change floating-point rounding. These
options are candidates for measurement, not a guarantee of faster solves or
bit-identical outputs. Check physical residuals, currents, and far fields against
the default build at the same numerical settings before adopting a candidate.
`-Ofast` and `-ffast-math` are not enabled by these switches because their
relaxed floating-point assumptions can affect solver arithmetic and checks.
See GCC's [optimization options](https://gcc.gnu.org/onlinedocs/gcc-11.4.0/gcc/Optimize-Options.html)
and [x86 options](https://gcc.gnu.org/onlinedocs/gcc-11.4.0/gcc/x86-Options.html).

Profile-guided optimization requires a separate training stage:

1. Build the instrumented executable:

   ```sh
   make -j2 PGO=generate BUILD_DIR=build/pgo TARGET=build/pgo/run_mlfma
   ```

2. Run that executable on representative inputs using the separate working
   directories described below. Profiles go to the absolute `PROFILE_DIR`, even
   when the solver runs elsewhere. Training is never performed by `make`.
   Agent-executed training requires explicit approval under `AGENTS.md`.
3. Preserve the instrumented executable and rebuild using the collected data:

   ```sh
   cp build/pgo/run_mlfma build/pgo/run_mlfma.instrumented
   make -j2 PGO=use BUILD_DIR=build/pgo TARGET=build/pgo/run_mlfma
   ```

Keep the same object directory, compiler, sources, and optimization settings
between collection and use; change only `PGO`. Use a fresh profile directory
when sources or settings change. Updated `.gcda` data triggers recompilation
in `PGO=use` mode, including nested profile directories from absolute object
paths. `make clean` retains profile data. Instrumented execution has
profiling overhead and should not be used as the optimized timing result.
See GCC's [instrumentation options](https://gcc.gnu.org/onlinedocs/gcc-11.4.0/gcc/Instrumentation-Options.html).

## Run a simulation

The executable takes no command-line arguments. It reads **input.pnu from the
current working directory** and writes results there. The mesh is embedded in
the input; `model_name` is a label, not an external mesh filename.

Use a fresh working directory to preserve earlier results. For the tiny
regression fixture (two octahedra, 24 RWGs / 48 unknowns, three illuminations):

```sh
repo_dir="$PWD"
run_dir=$(mktemp -d /tmp/mlfma-ibc-run.XXXXXX)
cp tests/ibc_reference.pnu "$run_dir/input.pnu"
(
  cd "$run_dir" || exit 1
  OMP_DYNAMIC=FALSE "$repo_dir/run_mlfma" > solver.log 2>&1
)
solver_status=$?
printf 'Exit status: %s\nResults: %s\n' "$solver_status" "$run_dir"
```

The original production inputs are retained unchanged:

| Input | Geometry | Original settings |
| --- | --- | --- |
| `src/input.pnu` | SOM, 91,923 points / 183,842 triangles | 5 GHz, mode 2, eight requested threads, principal cuts |
| `src/input_sphere.pnu` | Sphere, 12,531 points / 25,058 triangles | 3 GHz, mode 2, eight requested threads, monostatic field |

Read each header before running; these are substantially larger than the test
fixture. Change settings in copies, and inspect both the exit status and
`result_solver_status.txt` before interpreting results.

## Input format and settings

The text reader expects the named parameter records in the order below.
An optional preamble precedes `SimulationInputs`; blank lines are allowed.
Do not add inline comments to parameter or mesh records: their field counts
are checked exactly. Coordinates are in metres, frequency in hertz, and input
angles in degrees. Theta is the polar angle from +z; phi is the azimuth about z.

| Record | Meaning |
| --- | --- |
| `model_name name` | Single-token model label, at most 199 bytes. |
| `Formulations 1` | Coupled IBC formulation; only 1 is implemented here. |
| `Frequency start end count` | Positive frequency endpoints; count must be 1. The start frequency is used. |
| `source_type 0` | Plane-wave illumination; only 0 is supported here. |
| `source_PW theta_start theta_end n_theta phi_start phi_end n_phi pol` | Incident propagation directions sampled inclusively. Counts must be positive. A count of 1 uses the corresponding start angle. `pol` is 0 or 1. |
| `Preconditioners mode` | 0: GMRES without preconditioning; 1: GMRES with left block-diagonal preconditioning using pivoted LU factors; 2: FGMRES with approximate-MLFMA FGMRES, itself preconditioned by near-only GMRES. |
| `num_mlfmalevels levels` | 1 selects the tree depth automatically; larger values request an explicit depth within checked index limits. |
| `dir_ff selector` | Far-field output mode; see below. |
| `num_digits value` | Positive parameter controlling MLFMA truncation; it is not a guarantee of that many accurate output digits. |
| `val_approx value` | Mode-2 harmonic reduction factor in `[0,1]`. A value of 0 selects a near-only inner operator; positive reduced grids must support the interpolation stencil. |
| `num_interpol value` | Positive angular interpolation stencil width; must fit every grid where it is used. |
| `val_itermax_pc value` | Positive inner GMRES iteration limit for mode 2. The approximate FGMRES tolerance is `0.5`; its near-only GMRES preconditioner keeps a 20-iteration budget and tolerance `0.4`. |
| `val_soltol value` | Outer relative residual tolerance, strictly between 0 and 1. |
| `dir_cfc theta phi` | Direction for the post-solve current/field contribution diagnostic. |
| `val_itermax value` | Positive outer iteration budget; production FGMRES also uses it as the restart length. |
| `parallel_omp enabled threads` | `enabled` is 0 or 1; requested thread count must be positive even when disabled. |

The solver sets the requested OpenMP team size from `parallel_omp`; setting
`OMP_NUM_THREADS` alone does not replace that input. Runtime limits such as
`OMP_THREAD_LIMIT` and dynamic adjustment can reduce the actual team. Use
`OMP_DYNAMIC=FALSE` for repeatable team requests. The solver handles actual teams
that differ from the request, but speedup depends on the case and hardware.

For polarization, `calculate_dir_pw()` constructs the electric-field direction
as `(theta_E, phi_E) = (fmod(theta + 90, 180), phi)` for `pol = 0`, and
`(90, fmod(phi + 90, 360))` for `pol = 1`, using C++ `fmod`.
This states the implementation's sign
convention explicitly. Illuminations are ordered with theta outermost and phi
innermost; result indices start at zero.

The mode-1 preconditioner factors each leaf's A11 and A22 diagonal blocks
separately, using partial-pivoted complex LU. This preserves the original
uncoupled preconditioning method; it does not factor a combined 2-current block.
No LU factors are constructed for modes 0 or 2.

Mesh records follow the parameter section:

```text
## End-Param
ModelParts <number_of_parts>
<part_name_0> <impedance_real_0> <impedance_imag_0>
... one name and complex impedance per part ...
POINTS <number_of_points>
<x_0> <y_0> <z_0>
... one coordinate record per point ...
CELLS <number_of_triangles>
<part_index> <node_0> <node_1> <node_2>
... one connectivity record per triangle ...
```

Part/node indices are zero based. There must be 1–50 parts. The finite, nonzero
complex material values are used directly in `sqrt(zeta)` and `sqrt(1/zeta)`
terms in `calculate_Znearmatrix()`; no new impedance conversion is introduced.
Zero impedance is unsupported because the implemented equations contain its
reciprocal. Reciprocal/square-root overflow is also rejected.

The reader validates finite coordinates, indices, nonzero edges/areas, and
duplicate triangles. Connectivity checks reject nonmanifold edges/vertices.
Use consistently oriented triangles; validation does not repair winding or
establish the physical suitability of IBC on an open target. RWGs span shared
edges; boundary edges do not create half-RWG unknowns.

`dir_ff` keeps this repository's original direct-field interpretation:

| Value | Output |
| --- | --- |
| `< 0` | Three principal-plane cuts (ZX, ZY, XY), 361 directions each. |
| `0` | Monostatic backscatter for each illumination. |
| `> 0` | Direct angular grid with this step in degrees, clamped to at least 0.1. The effective step must divide 180 and 360; incompatible allocations are rejected. |

There is no accelerated far-field grid selector in this IBC executable.

## Outputs and convergence

Files retain the historical spelling **farefield**.

| File | Contents |
| --- | --- |
| `result_solver_status.txt` | Per-illumination overall/Krylov status, iterations, outer product count, true Krylov and physical residuals, tolerance. |
| `result_farefield.txt` | Monostatic results, appended in illumination order within the run. |
| `result_farefield_<index>.txt` | Direct cuts or angular grid for each illumination. |
| `result_rwgcoefficient.txt` | Final illumination's 2N complex coefficients: N scaled electric, then N magnetic. |
| `result_currentontris_J.txt` | Final illumination's electric-current magnitude per triangle. |
| `result_currentontris_M.txt` | Final illumination's magnetic-current magnitude per triangle. |
| `result_cfcontris.txt` | Leaf-box/total field contribution magnitude in `dir_cfc`, assigned to triangles. |

Far-field columns are theta, phi, `abs(E_theta)`, `phase(E_theta)`,
`abs(E_phi)`, `phase(E_phi)`: directions in degrees, phases in radians.
They are far-zone field coefficients with radial propagation omitted,
not RCS in dB. Normalization remains in `calculate_farfield()`.

The physical residual is `||b - A*x|| / ||b||`, evaluated over **both currents**
using the original RHS and accurate, unpreconditioned MLFMA operator. A zero RHS
uses an absolute residual. Mode 1 also reports the residual of the
left-preconditioned Krylov equation. Its convergence does not ensure that the
physical tolerance is met: the tiny regression fixture intentionally checks this
failure status for mode 1. Outer product counts exclude inner solves and the
additional physical-residual check.

Exit 0 means every illumination converged and the status file opened; exit 1
indicates an input error or reported solve/status-file failure. Other output
files may still be written for failed solves. Residual convergence verifies
the implemented discrete equation, not mesh/quadrature/MLFMA accuracy.

## Code layout and solution path

1. `read_simulationinput()` validates parameters/materials and the embedded mesh.
2. `generate_edgeconnection()` builds RWGs and triangle geometry.
3. `generate_octtree()`, `generate_nearfarboxes()`, `generate_Z*()`, and
   `generate_translationcomponent()` allocate interaction and hierarchy data.
4. `calculate_loadbalance()`, `calculate_Znearmatrix()`, and
   `calculate_translationoperator()` prepare partitions, four near blocks, and
   far translations; mode 1 additionally factors its two diagonal block sets.
5. `stage_solution()` builds each 2N RHS, runs GMRES/FGMRES, checks the physical
   residual, and writes each illumination's far field.
6. `stage_aftersolution()` writes final currents/diagnostics; `clear_memory()`
   releases all owned simulation data and cached workspace.

`mvm.cpp` combines the four-block near product (`mvm_Zneardotx.cpp`) with
`mvm_Zfar1121dotx.cpp` and `mvm_Zfar1222dotx.cpp`. Their `_in` variants provide
approximate products for mode 2. `krylov.cpp`, `gmres.cpp`, and `fgmres*.cpp`
implement the solvers. `mlfma_workspace.*` owns reusable worker scratch and
shared angular/stencil/phase tables. Global state makes this a standalone,
non-reentrant solver.

## Validation and development history

See [tests/README.md](tests/README.md) for the bounded IBC checks and
[PORT_NOTES.md](PORT_NOTES.md) for adaptations from the PEC repository.
Builds never execute tests. Agent-executed simulations/tests/benchmarks require
the approval specified in AGENTS.md.

[CHANGELOG_DEV.md](CHANGELOG_DEV.md) records changes and checks;
[HANDOFF.md](HANDOFF.md) records current status and artifact locations.

## License

Project-owned code and documentation use the [MIT License](LICENSE).
`src/times.hpp` retains Mehmet Nuri's existing copyright notice; the project
license does not override third-party rights for that separately attributed file.
