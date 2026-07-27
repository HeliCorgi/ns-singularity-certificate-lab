# Reproducibility protocol

## Minimal environment

- CPython 3.10 or newer;
- NumPy;
- pytest for validation;
- no GPU and no network access at run time.

Create an isolated environment if desired, then install the repository in
editable mode:

```console
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

On POSIX shells, use `source .venv/bin/activate` before the two `python`
commands.  The replay below assumes that this environment remains active.

## Required replay sequence

```console
python -m pytest
python experiments/run_manufactured.py ^
  --config configs/manufactured.json ^
  --output-dir outputs/manufactured_replay
python experiments/run_baseline.py ^
  --config configs/baseline.json ^
  --output-dir outputs/baseline_replay
python -m experiments.run_time_convergence ^
  --config configs/baseline_time_convergence.json ^
  --output-dir outputs/time_convergence_replay
```

The line continuations above are for Windows `cmd.exe`; PowerShell users may
place each command on one line.  Use new output directory names so existing
evidence is not overwritten.

## What to compare

1. Every test passes, including the intentional fault-injection suite.
2. Manufactured diagnostics show decreasing error and the documented observed
   orders for all six measured identities, including the independent
   Cartesian-reconstruction/curl route.
3. The test suite's uniform-Cartesian manufactured field shows decreasing
   approximately second-order errors for three-component divergence, full
   curl, vector Laplacian, primitive residual terms, and total residual.  Its
   saved-candidate test performs a real artifact round-trip and rejects each
   intentional cylindrical fault.
4. Baseline `accepted_as_control` is true, energy is nonincreasing, the
   boundary comparison passes, and the false-positive guard refuses a blow-up
   fit.
5. Time-convergence `accepted_as_temporal_control` is true; the actual steps
   are exactly \(0.5,0.25,0.125\) on one identical \(513\)-point grid, analytic
   errors decrease, both raw-error orders exceed 2, the step-doubling order is
   approximately 2, energy does not increase, and the boundary test passes.
6. Candidate, run-config, and diagnostic loaders accept the replay artifacts;
   all three contain an identical pre-write provenance snapshot.
7. Verify `manifest.json.sha256` for the time-convergence output and the
   per-file SHA-256 values inside its manifest before reading the data as
   evidence.
8. Compare physical numeric values with tolerances, not compressed NPZ bytes:
   ZIP metadata and floating-point library details can change bytewise hashes
   across platforms.

The v2 writers record Python, NumPy, operating system, Git state, and a stable
source-input fingerprint automatically.  Record any environmental deviation
that those fields do not capture.  Current reference directories are
`outputs/manufactured_v5`, `outputs/baseline_v5`, and
`outputs/time_convergence_v1`; replay into new names because all experiment
scripts refuse to overwrite non-empty evidence.  A numerically close replay
is not an independent implementation and is not a mathematical proof; a
future candidate requires separately written discretizations, rigorous error
bounds, and further obligations listed in `docs/proof_obligations.md`.
