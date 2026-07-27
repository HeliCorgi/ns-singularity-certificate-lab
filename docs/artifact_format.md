# Artifact and diagnostic formats

## Explicit candidate bundle

`save_candidate("name.npz", ...)` creates three files:

1. `name.npz` — explicit canonical little-endian binary64 (`<f8`) arrays
   `r`, `z`, and named fields;
2. `name.manifest.json` — schema version, array shapes, dtypes, per-array
   SHA-256 values, archive SHA-256, grid metadata, configuration, seed, and
   descriptive metadata;
3. `name.manifest.json.sha256` — integrity value for the manifest itself.

The current schema identifier is `ns-certificate-lab/candidate/v2`.  Its
required `description` object records:

- representation and coordinate system;
- units for `r`, `z`, and every stored field;
- nondimensionalization/normalization;
- physical time and positive viscosity;
- basis, array-order, velocity-recovery, and elliptic-sign convention.

The required `provenance` object records the Python implementation/version,
NumPy version, platform, Git HEAD or an explicit unborn/not-a-checkout
sentinel, dirty state, and a stable SHA-256 fingerprint over
`pyproject.toml`, `src/`, `experiments/`, and `configs/`.  Generated
`outputs/`, caches, bytecode, and editable-install `.egg-info` metadata are
excluded from the source fingerprint.  A logical experiment takes this
snapshot once before writing and injects the same object into its run
configuration, diagnostics, and candidate manifest.

Loading uses `allow_pickle=False` and rejects missing files, nonfinite arrays,
non-`<f8` arrays (even when a manifest has been recomputed), unknown schemas,
malformed names, inconsistent shapes, invalid provenance, and any checksum
mismatch.  JSON parsing also rejects duplicate keys, `NaN`, infinities, and
finite-parser overflow.

Legacy `ns-certificate-lab/candidate/v1` bundles remain read-only compatible
when their arrays satisfy the same canonical dtype and integrity checks.  They
load with `description=None`, `provenance=None`, and
`legacy_provenance_absent=True`; the absence is never silently filled.

A neural-network state is intentionally not a valid candidate field.  A
network-assisted search must first evaluate and export the retained object as
explicit arrays or, in a future schema, explicit basis coefficients.

## Run configuration

`save_run_config` writes `ns-certificate-lab/run-config/v2` canonical JSON
containing a concrete integer seed, finite JSON parameters, and the runtime
provenance snapshot, plus a SHA-256 sidecar.  Version 1 is accepted only as a
legacy read with the same explicit provenance-absence flag.
`seed_everything` seeds both Python and NumPy and returns an explicit NumPy
generator.

## Diagnostics

The JSON diagnostic schema is
`ns-certificate-lab/diagnostics/v1`; CSV diagnostics use scalar columns.  Both
formats receive checksum sidecars, and their readers verify integrity before
parsing.  Manufactured JSON additionally stores the full axis-regularity
report and RMS values for every signed PDE term, forced residual error, and
independent Cartesian-reconstruction/curl defect.  The test suite deliberately
modifies a recorded value and requires the read to fail.

The fixed-grid temporal control uses a directory-level evidence bundle rather
than the candidate schema.  It contains the exact config snapshot,
`summary.json`, per-step and per-run CSV files, and explicit final profiles in
an `allow_pickle=False` NPZ.  `manifest.json` fixes the expected five payload
files by SHA-256 and byte length, while `manifest.json.sha256` protects the
manifest.  Its auxiliary \(R=3,4\) boundary-radius comparison is labelled
separately because it is not a direct truncation-error bound for the main
\(R=5\) run.

## Security and proof limitations

SHA-256 detects accidental or deliberate mutation relative to the recorded
manifest.  It does not establish who created the artifact, that the program was
correct, or that floating-point values enclose exact real quantities.  Future
proof certificates need signed provenance where appropriate and interval or
exact-arithmetic enclosures for every proof-relevant bound.
