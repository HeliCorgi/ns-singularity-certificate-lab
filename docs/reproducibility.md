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
```

The line continuations above are for Windows `cmd.exe`; PowerShell users may
place each command on one line.  Use new output directory names so existing
evidence is not overwritten.

## What to compare

1. Every test passes, including the intentional fault-injection suite.
2. Manufactured diagnostics show decreasing error and the documented observed
   orders for all six measured identities, including the independent
   Cartesian-reconstruction/curl route.
3. Baseline `accepted_as_control` is true, energy is nonincreasing, the
   boundary comparison passes, and the false-positive guard refuses a blow-up
   fit.
4. Candidate, run-config, and diagnostic loaders accept the replay artifacts;
   all three contain an identical pre-write provenance snapshot.
5. Compare physical numeric values with tolerances, not compressed NPZ bytes:
   ZIP metadata and floating-point library details can change bytewise hashes
   across platforms.

The v2 writers record Python, NumPy, operating system, Git state, and a stable
source-input fingerprint automatically.  Record any environmental deviation
that those fields do not capture.  A numerically close replay is not an
independent implementation; a future candidate requires a separately written
discretization and diagnostics.
