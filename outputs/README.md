# Generated outputs

Small, reviewable outputs from the manufactured audit, non-singular spatial
control, and fixed-grid temporal control are tracked here.  Every output set
must contain its exact input configuration and a machine-readable summary.
Candidate-search outputs must use a distinct run identifier and must not
overwrite prior evidence.

- `manufactured/` and `baseline/` are the preserved initial artifacts.
  `manufactured/` uses the legacy candidate/run-config v1 format and therefore
  has no runtime provenance.
- `manufactured_v2/`, `baseline_v2/`, `manufactured_v3/`, and `baseline_v3/`
  are preserved intermediate provenance-format runs made before all
  normalization and editable-install fingerprint fixes were complete.
- `manufactured_v4/` and `baseline_v4/` are preserved evidence from before the
  uniform-Cartesian and fixed-time-step audit code was added.
- `manufactured_v5/` is the current manufactured audit replayed against the
  present source fingerprint.  It has candidate v2, shared source/runtime
  provenance, full transformed-PDE term diagnostics, axis metrics, and the
  independent cylindrical-sampling Cartesian-reconstruction/curl check.
- `baseline_v5/` is the current non-singular spatial-refinement control with
  explicit normalization and units, a provenance-bearing summary, and an
  integrity sidecar for its file manifest.
- `time_convergence_v1/` is the current fixed-\(n_r=513\) temporal control at
  \(\Delta t,\Delta t/2,\Delta t/4\).  `summary.json` and
  `time_convergence.csv` record analytic and step-doubling errors/orders,
  energy, maximum physical vorticity, and boundary sensitivity;
  `time_diagnostics.csv` stores histories, `final_profiles.npz` stores
  explicit final arrays, and the manifest has a digest sidecar.

The uniform-Cartesian divergence/full-curl/primitive-residual route and the
saved-candidate end-to-end fault injections are deterministic automated tests;
they do not currently create a separate evidence directory.

The old directories are retained deliberately; adverse or superseded evidence
must not be silently rewritten or deleted.
