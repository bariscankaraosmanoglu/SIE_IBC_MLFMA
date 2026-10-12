# Bounded IBC regression checks

`check_port.py` compiles standalone checks and runs only tiny fixtures, retaining
all artifacts in the specified directory (default: unique `/tmp` directory).
Python 3, g++, and OpenMP suffice; no external numerical packages are needed.
An agent must obtain AGENTS.md approval before runtime execution.

```sh
python3 -B tests/check_port.py --output-dir /tmp/<unique-approved-directory>
python3 -B tests/check_port.py --sanitize --output-dir /tmp/<other-approved-directory>
```

`--build-only` compiles without running anything. `--skip-build` uses binaries
already compiled from the same sources. `--baseline-source /path/to/original/src`
optionally checks guarded original operators in a normal build; it writes its
guarded copy only to the output directory, never edits the original source.

The suite covers complex Krylov rotations/breakdown/statuses, preconditioner
failure propagation, physical-residual reporting, pivoted LU and aliasing,
angular/stencil/phase reuse, surface topology, and thread partitions. Coupled
IBC checks use `ibc_reference.pnu`: two closed octahedra, 12 points, 16 triangles,
24 RWGs / 48 unknowns; 1 GHz; complex materials 0.001+0i and 0.2+0.1i;
three PW illuminations (theta 37, phi 19–21); digits 6; interpolation 4;
approximation 0.8; outer tolerance 0.001, budget 100; approximate FGMRES budget 10 / tolerance 0.5; near-only GMRES budget 20 / tolerance 0.4.

Actual operator cases request 1, 2, 8, 16, and 32 workers, plus disabled,
dynamic, and limited teams. Tests compare all four near matrices, both far
branches, and approximate products to serial at relative tolerance 1e-12;
12 repeated products check stable workspace capacity. A test-only fully assembled
direct-MoM reference checks MLFMA accuracy at 0.001. A test-only dense accurate
48-by-48 MLFMA operator independently checks final 2N solver residuals for
modes 0/1/2. Mode 1 on this fixture converges only in the preconditioned
equation: the physical residual exceeds 0.001. The suite requires the correct
`physical_residual_not_converged` status and exit 1, independently verifies the
residual, and keeps the tolerance unchanged. Modes 0/2 must converge physically. Three illuminations exercise RHS reuse/ownership. A zero-approximation
case checks near-only inner solves. Sixteen malformed inputs must exit 1 with
an input diagnostic. ASan/UBSan repeat the suite; a final operator case checks
leaks and repeated cleanup/rebuild.

The guarded original reference fixes only the absent-basis sentinel lookup and
uses safe serial thread ranges; original numerical kernels are otherwise
unchanged. The suite does not execute either production input, measure scaling,
or establish analytical IBC accuracy. Logs, commands.json, fixture copies,
binary dumps, and results.json are retained. Failed attempts must also be kept.
