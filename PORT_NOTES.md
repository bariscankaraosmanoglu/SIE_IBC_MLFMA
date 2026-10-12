# PEC-to-IBC port

Source: SIE_PEC_MLFMA `d52f36c886ba1a7a5ad93688c72db4cf773808a7`.
Target starting point: SIE_IBC_MLFMA `d3215e11baaffd1b6368bed473261b14759ad11d`,
branch `codex/code_adaptation`. No pre-existing target working-tree changes.
The source repository and both original IBC production inputs remain unchanged.

## Changes carried over

| Prior item | IBC implementation |
| --- | --- |
| 1: thread partitions | `calculate_loadbalance.*`, near assembly/product, and all four far kernels use actual OpenMP teams; weighted near prefixes replace fixed worker arrays. Handles disabled, limited, dynamic, and >16-thread teams. |
| 2: coefficient/declaration correctness | `solver_globals.h` gives shared arrays/index pointers consistent types. PEC CFIE coefficient forcing is intentionally inapplicable to IBC. Connectivity validation is shared without changing IBC material terms. |
| 3: Krylov correctness | `krylov.*`, `gmres.cpp`, `fgmres*.cpp`, and preconditioner wrappers handle complex rotations, bounds, breakdown, and termination. `stage_solution()` uses 2N dimensions and reports the original coupled physical residual. |
| 4: input checks | `read_simulationinput.cpp`, `input_validation.*`, allocation/tree modules reject malformed and unsupported inputs before computation. Retain all parameter names, complex material records, and direct IBC `dir_ff` semantics. |
| 5: block factors | `matrixBlock.*`, `calculate_blockdiagonalpreconditioner.cpp`, `mvm_BDPmatrix.cpp` replace explicit inverses with pivoted LU, mode 1 only. Preserve separate A11/A22 application; remove unused A12/A21 inverse construction. |
| 6: ownership/workspace | `clear_memory.cpp`, `mlfma_workspace.*`, interpolation and four far kernels own/free data consistently, reuse worker scratch, and allow repeated cleanup/rebuild. Krylov storage redesign remains excluded. |
| 7: shared caches | Shared angular grids, interpolation stencils, and eight child shifts per level serve accurate and positive approximate far operators. No full per-RWG radiation cache. |
| Build/docs/license | Incremental Makefile with unchanged default flags, optional LTO/CPU/PGO switches, IBC README, MIT license and author notice, development notes, bounded tests. |

## Confirmed IBC-specific defects addressed

- `calculate_Znearmatrix()` indexed `tri_media[basistri]` before checking the
  `-1` sentinel for an absent neighboring basis. Move the material lookup inside
  the existing guard; retain all four equation terms and quadrature.
- The original fixed thread partition builder writes beyond its radiation
  range allocation. Its replacement uses checked ranges from the actual team.
- Global declarations disagreed on array shapes/index types; centralize them.
- Original block setup formed four inverse sets although application used only
  A11/A22. Factor only those sets, and only when mode 1 is selected.
- Original cleanup left several IBC allocations owned or dangling. Release all
  four near matrices and both current block factors, approximate translations,
  and workspace; reset owned pointers.

The IBC operator has 2N unknowns. Mode 2 retains its original nested structure:
outer accurate FGMRES → approximate FGMRES (input inner budget, tolerance 0.5) →
near-only GMRES (20 iterations, tolerance 0.4). This differs from the PEC
preconditioner and must not be replaced by its single GMRES/0.1 policy.
All inner solution entries are initialized across both current blocks.

The existing signs, source/receiving kernels,
`(i*k/(4*pi))^2` prefactor, material square roots, current scaling, and direct
field formulas are retained. Applying PEC open-surface EFIE weights here would
change the intended equations and was not done.

- `calculate_translationoperator()` initialized P1 even when its expansion
  contained only P0. Guard that initialization for zero-order expansions. The
  approved runtime cases exercise positive orders and zero-approximation
  near-only solves, not the positive-factor/zero-order translation branch.

## Evidence and limits

The test driver compares near/far/approximate products with original IBC
kernels using two explicitly documented guards: safe material-sentinel lookup
and a serial partition allocation. This is a **guarded original operator
reference**, not a claim that the unmodified old executable runs correctly.
The original source and unmodified compiled executable are retained separately.

Small tests materialize a 48-by-48 operator only for independent residual
checks and compare against a fully assembled direct-MoM test matrix. No dense
production-size matrix is introduced. See HANDOFF.md for executed results.

No analytical IBC scattering reference, production SOM/sphere simulation,
large-core performance benchmark, or optional-compiler-option benchmark is
included in the bounded plan. Performance benefits are expected, not measured.
Open-surface IBC physical suitability and winding/self-intersections remain
outside these checks.

Final bounded normal and ASan/UBSan suites passed (45 cases each, with the
sanitizer suite including its leak check). Thread configurations gave identical
products. The largest guarded-original relative difference was 3.34e-15; the
small direct-MoM/MLFMA comparison error was 3.91e-9 at tolerance 0.001.
Modes 0 and 2 converged physically. Mode 1 correctly reported physical residuals
0.00120–0.00168 despite preconditioned convergence; the tolerance was unchanged.

Performance and production-input validation remain unmeasured. A dedicated
positive-factor/zero-order translation run and a two-level-tree run were not
part of the approved settings; the additional allocation/bounds guards were
reviewed and compiled, with the normal approved configurations rechecked.
