# Development change log

### 2026-10-12 21:00 UTC — Port PEC infrastructure improvements to coupled IBC
- User request and authorized scope: analyze IBC variant and apply applicable PEC
  improvements; source/docs/build/test edits authorized in this repository.
- Branch and starting HEAD: codex/code_adaptation, d3215e1; PEC source d52f36c.
- Files and functions changed: see PORT_NOTES.md mapping and the review diff;
  actual-team partitions, centralized declarations, Krylov and stage_solution,
  input/allocation validation, two-current pivoted block factors, ownership,
  shared MLFMA workspace/caches, Makefile, README, tests, license.
- Changes and rationale: preserve the coupled 2N IBC operator; fix confirmed
  absent-basis material lookup and thread-array bounds defects. Reuse applicable
  PEC correctness/performance work without applying PEC CFIE/EFIE rules.
- Expected numerical impact: operator terms unchanged except removal of undefined
  accesses; corrected Krylov and pivoted factors can change convergence/rounding.
- Validation plan: bounded normal and ASan/UBSan check_port.py suites, two
  octahedra/48 unknowns/1 GHz/three illuminations/tolerance .001, modes 0/1/2,
  actual 1–32-thread products and malformed inputs; guarded original comparison.
- Approval received and its scope: explicit user “approved” for both bounded
  suites; no production sphere/SOM simulation or performance measurement.
- Commands actually executed: original g++ full build; prepared Makefile -j2;
  check_port.py --build-only --output-dir /tmp/mlfma-ibc-port.373jot77/validation
  --baseline-source /tmp/mlfma-ibc-port.373jot77/before/src; detailed build commands
  in that directory's commands.json. Read-only source/Git inspection.
- Build outcome: original full build and prepared default Makefile/test builds
  passed; initial missing val_approx declaration corrected, failed log retained.
- Test outcome: NOT RUN yet — explicit approval now received.
- Measurements and output locations: /tmp/mlfma-ibc-port.373jot77; no runtime or
  RSS comparison measured; all before/prepared/validation artifacts retained.
- Remaining uncertainty or pending work: run approved checks, review final diff;
  production performance and analytical IBC accuracy remain unmeasured.
- Commit/push status: neither requested; changes remain uncommitted/unpushed.

Final implementation and validation update — 2026-10-12 21:46 UTC:
- IBC review adaptations: preserve its original mode-2 nested approximate FGMRES
  (input budget, tolerance 0.5) with near-only GMRES (20 iterations, tolerance
  0.4); initialize both current blocks and propagate inner failures. Preserve
  compact active-level aggregation storage, including a safe two-level leaf.
  Correct remaining RHS array declarations/2N size diagnostic. Guard P1
  initialization when only P0 is allocated in translation setup.
- Commands actually executed (repository cwd /home/bariscank/repo/SIE_IBC_MLFMA):
  `make -j2 BUILD_DIR=/tmp/mlfma-ibc-port.373jot77/target-build
  TARGET=/tmp/mlfma-ibc-port.373jot77/run_mlfma`; repeated incrementally after
  focused fixes. Normal and sanitizer check_port.py build/runs in validation,
  validation-san, validation-final, validation-san-final, and reviewed directories;
  all detailed compiler/runtime commands and exit codes retained in commands.json.
  Final normal command:
  `python3 -B tests/check_port.py --skip-build --output-dir
  /tmp/mlfma-ibc-port.373jot77/validation-reviewed --baseline-source
  /tmp/mlfma-ibc-port.373jot77/before/src`.
  Final sanitizer build:
  `make -j2 BUILD_DIR=/tmp/mlfma-ibc-port.373jot77/reviewed-san-build
  TARGET=/tmp/mlfma-ibc-port.373jot77/reviewed-san-run_mlfma
  CXXFLAGS='-Wall -fexceptions -O1 -g -std=gnu++17
  -fsanitize=address,undefined -fno-omit-frame-pointer'`.
  `python3 -B /tmp/mlfma-ibc-port.373jot77/build_reviewed_drivers.py --sanitize`
  links the instrumented operator driver with those objects; unchanged unit
  binaries/provenance are retained from the final suite builds.
  Final sanitizer command:
  `python3 -B tests/check_port.py --sanitize --skip-build --output-dir
  /tmp/mlfma-ibc-port.373jot77/validation-san-reviewed`.
- Build outcome: default O3/GNU++17/OpenMP and O1 ASan/UBSan builds pass with
  GCC 11.4.0 / GNU Make 4.3. Existing unused-variable warnings remain.
- Test outcome: final normal 45/45; final sanitizer 45/45 including LeakSanitizer.
  All thread products equal serial exactly; guarded-original max relative
  difference 3.337320896378312e-15. Direct-MoM/MLFMA relative difference
  3.9135508153203795e-9 <= unchanged .001 tolerance. Pivoted block application
  residual 2.138244283071602e-16. Repeated workspace capacity/cleanup checks pass.
  Modes 0 and 2 physically converge; mode 1 correctly exits 1 with
  physical_residual_not_converged (0.0012033–0.0016821 at tolerance .001).
  Independent 48-by-48 operator checks confirm each final reported residual
  within 1e-12. All 16 invalid inputs reject as expected; near-only inner solves
  converge. Original IBC mode-2 budgets/tolerances were retained.
- Failed/intermediate attempts retained: omitted block-test argument; initial
  test assumption that left-preconditioned convergence guaranteed physical
  convergence; superseded PEC-style inner preconditioning; LeakSanitizer's
  sandbox ptrace failure. Final leak check passed outside sandbox tracing using
  the same bounded approved cases. No actual leak/bounds/UB report remains.
- Output locations: /tmp/mlfma-ibc-port.373jot77, especially validation-reviewed,
  validation-san-reviewed, validation-summary.json, before, prepared, *build*.log,
  normal-*.log, sanitizer-*.log, and archived intermediate attempts.
- Measurements: no before/after runtime or peak-RSS benchmark performed; no
  production performance benefit claimed. Stable workspace capacity is a test
  check, not an RSS measurement. Original large inputs retained byte for byte.
- NOT RUN: production SOM/sphere cases, analytical IBC reference, large-core
  scaling, optional compiler-flag benchmarks, dedicated positive-factor/zero-order
  translation and two-level-tree boundary cases. Expanded runtime cases require
  new AGENTS.md approval. Static guards for the last two cases are compiled but
  their specific branches were not runtime-validated in this plan.
- Final review: applicable physics/signs/current scaling and nested settings
  retained; source PEC clean; target production inputs/times.hpp/AGENTS.md
  unchanged; tracked diff whitespace check passes with existing CRLF handling.
  New source, test, build and documentation files reviewed separately.
- Commit/push status: all changes remain uncommitted and unpushed on
  codex/code_adaptation, HEAD d3215e1. No commit/push authorization received.


### 2026-10-12 21:51 UTC — Prepare reviewed IBC port commit
- User request and authorized scope: commit the completed IBC port on the current
  branch; user will push manually. No push requested or performed.
- Branch and starting HEAD: codex/code_adaptation, d3215e1.
- Files and functions changed: development notes updated for commit authorization;
  source/build/test files match the final validated source manifest.
- Changes and rationale: package the reviewed 64-file port with its validation
  evidence, preserving both original production inputs and all retained artifacts.
- Expected numerical impact: none from this commit preparation; documentation only.
- Validation plan: verify manifest, selected staged paths and staged whitespace;
  reuse completed normal 45/45 and ASan/UBSan 45/45 checks. No tests rerun.
- Approval received and its scope: explicit user request to commit, excluding push.
- Commands actually executed: read-only Git status/branch/HEAD/config/hooks review;
  source SHA-256 comparison against final-source-manifest.json. Commit preparation
  stages only the 64 reviewed port paths and inspects the cached diff.
- Build outcome: unchanged, default and instrumented builds previously passed.
- Test outcome: unchanged, previously recorded bounded suites passed; known mode-1
  physical-residual failures and untested boundary/production cases remain recorded.
- Measurements and output locations: /tmp/mlfma-ibc-port.373jot77 retained;
  commit-review artifacts in /tmp/mlfma-ibc-commit.pwm76S.
- Remaining uncertainty or pending work: production performance/analytical accuracy
  and dedicated boundary runs remain outside the approved validation scope.
- Commit/push status: commit authorized and prepared with message
  'Port solver correctness and workspace improvements to IBC'; this entry and
  handoff are included in that commit. The resulting identifier is reported to
  the user after creation. No push performed.
