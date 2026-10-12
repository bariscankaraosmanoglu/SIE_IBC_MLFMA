# Current handoff

Task complete: applicable PEC infrastructure improvements ported to coupled IBC,
with IBC-specific signs, complex material terms, 2N current ordering, direct
far-field modes, and original nested mode-2 preconditioner retained. See
PORT_NOTES.md for the mapping and confirmed fixes.

Branch: codex/code_adaptation. Starting HEAD d3215e1; PEC source d52f36c.
No pre-existing user changes. The user requested the reviewed port commit; this
handoff is included in that commit. No push is performed; the user will push.
The resulting commit identifier is reported after creation. Source PEC, original target production inputs, times.hpp
notice, and AGENTS.md remain unchanged. All artifacts retained.

## Completed checks

- Default Makefile and ASan/UBSan builds pass with GCC 11.4.0 / Make 4.3.
- Final normal suite 45/45; sanitizer suite 45/45 including leak detection.
- Serial/disabled/2/8/16/32/dynamic/limited-team coupled products agree exactly.
- Guarded-original operator max relative difference 3.34e-15. The original
  reference guards only the material sentinel and unsafe serial partition
  allocation; it is not an unmodified original-solver comparison.
- Direct-MoM/MLFMA difference 3.91e-9 at unchanged tolerance .001. Pivoted block
  application residual 2.14e-16; aliasing, mode-only construction, repeated
  cleanup/rebuild and stable workspace checks pass.
- Modes 0/2 physically converge for all three illuminations. Mode 1 correctly
  reports physical_residual_not_converged and exits 1, despite convergence of
  the preconditioned equation. Physical residuals are independently checked
  with a test-only 48-by-48 accurate operator; production stays MLFMA.
- Sixteen malformed inputs reject; near-only inner solve passes. Mode 2 keeps
  approximate FGMRES tolerance .5/input budget, then near-only GMRES tolerance
  .4/budget 20. Both current blocks are initialized; inner failures propagate.
- Final diff and new files reviewed; whitespace check passes for existing CRLF.
  Initial argument omission, superseded preconditioner draft, and all failed
  sandbox LeakSanitizer attempts are retained and documented.

## Commands and artifacts

Build: `make -j2`. Actual build products were redirected to the retained task
folder with BUILD_DIR/TARGET overrides; no executable was left in the repo root.
Use `make -j2` locally to produce ./run_mlfma.

Final approved runtime commands (already executed):

```sh
python3 -B tests/check_port.py --skip-build --output-dir /tmp/mlfma-ibc-port.373jot77/validation-reviewed --baseline-source /tmp/mlfma-ibc-port.373jot77/before/src
python3 -B tests/check_port.py --sanitize --skip-build --output-dir /tmp/mlfma-ibc-port.373jot77/validation-san-reviewed
```

Artifacts: `/tmp/mlfma-ibc-port.373jot77`:
- before/: original sources and unmodified baseline executable (compiled only).
- prepared/: final staged source/docs/tests plus an earlier prepared executable;
  that earlier executable is superseded. Do not rerun preparation generators to
  overwrite final files; use the source manifest and review diff for provenance.
- target-build/, run_mlfma: final default build and executable.
- reviewed-san-build/, reviewed-san-run_mlfma: final instrumented build.
- validation-reviewed/, validation-san-reviewed/: final passing outputs, commands,
  result summaries and driver build provenance. The latter ran outside tracing
  because LeakSanitizer cannot operate under the sandbox's ptrace environment.
- validation-summary.json, metadata.json, applied_manifest.json: provenance.
- All preceding validation directories, snapshots and logs: retained attempts.

## Limits and next steps

Commit authorized by the user and prepared with the reviewed 64-file port.
The user will push manually. No agent push authorized. No full radiation-pattern
cache or separate Krylov storage redesign introduced. No production runtime/RSS
comparison measured; performance benefits are expected, not measured.

Production sphere/SOM, analytical IBC accuracy, optional compiler options, and
scaling beyond correctness checks at 32 workers were NOT RUN. Dedicated two-level
allocation and positive-factor/zero-order translation cases were NOT RUN; those
static guards are reviewed/compiled but their specific boundary branches need
an approved runtime plan. Do not change tolerances to hide mode-1 physical
residual failures. Improving that stopping criterion would be a separate task.
