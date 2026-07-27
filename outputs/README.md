# Generated outputs

Small, reviewable outputs from the manufactured audit and non-singular control
experiment are tracked here.  Every output set must contain its exact input
configuration and a machine-readable summary.  Candidate-search outputs must
use a distinct run identifier and must not overwrite prior evidence.

- `manufactured/` and `baseline/` are the preserved initial artifacts.
  `manufactured/` uses the legacy candidate/run-config v1 format and therefore
  has no runtime provenance.
- `manufactured_v2/`, `baseline_v2/`, `manufactured_v3/`, and `baseline_v3/`
  are preserved intermediate provenance-format runs made before all
  normalization and editable-install fingerprint fixes were complete.
- `manufactured_v4/` is the current manufactured audit with candidate v2,
  shared source/runtime provenance, full PDE-term diagnostics, axis metrics,
  and the independent Cartesian-reconstruction/curl check.
- `baseline_v4/` is the current non-singular control with explicit
  normalization and units, a provenance-bearing summary, and an integrity
  sidecar for its file manifest.

The old directories are retained deliberately; adverse or superseded evidence
must not be silently rewritten or deleted.
