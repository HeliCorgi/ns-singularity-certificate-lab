# Repository working agreement

This repository is a research scaffold for **candidate discovery and future
computer-assisted verification**, not a claim that a Navier--Stokes
singularity exists.

## Epistemic rules

1. Never label numerical growth, a fitted blow-up time, a small residual, or a
   neural-network output as a proof of singularity.
2. Trace every implemented reduced equation back to the three-dimensional
   incompressible Navier--Stokes equations.  Keep the physical three-dimensional
   divergence constraint separate from the formal five-dimensional scalar
   operator used by the regularized axisymmetric variables.
3. A formula marked `未確認`, `不整合`, or `誤り` in
   `docs/equation_audit.md` must not be an implementation premise.
4. Preserve negative results and failed runs.  New generated artifacts go in a
   new output directory or use a new run identifier; do not silently replace a
   preregistered result.
5. Prefer independent identities, manufactured solutions, convergence studies,
   checksums, and cross-implementation checks over self-consistency alone.
6. Record seeds, configuration, software versions, source revision, and
   checksums needed to replay every candidate or experiment.
7. Expensive searches require a separate, reviewed experiment specification.
   The default commands in this repository must remain CPU-friendly.

## Development rules

- Supported Python: 3.10 or newer.
- Install with `python -m pip install -e ".[dev]"`.
- Run the full test suite with `python -m pytest`.
- Run the non-singular benchmark with the command in `README.md`, always
  choosing a fresh `--output-dir`.
- Keep the numerical core in `src/ns_certificate_lab/`, tests in `tests/`,
  immutable experiment inputs in `configs/`, and generated evidence in
  `outputs/`.
- Candidate arrays and their manifests belong in `certificates/`; a manifest
  without matching hashes is invalid.
- Tests that intentionally inject a defect must assert that the defect is
  detected.  They must not be weakened merely to make CI green.
- If a derivation, theorem hypothesis, or citation cannot be checked, record it
  explicitly as unverified instead of filling the gap by inference.

## Review gates

An exploratory candidate may advance only after passing, in order:

1. equation audit;
2. axis regularity and boundary audit;
3. residual decomposition and refinement study;
4. independent reload and checksum verification;
5. domain/truncation and precision studies;
6. an independently implemented reproduction;
7. a separately scoped interval-arithmetic validation plan.

Passing these gates still produces a numerical candidate, not a counterexample.
