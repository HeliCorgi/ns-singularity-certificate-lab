# Reproducibility and release process

This document defines a release-grade workflow for the research scaffold.  A
release freezes an auditable numerical foundation; it is not evidence that a
Navier--Stokes singularity exists.

## Reproducibility levels

Three different claims must not be conflated:

1. **Artifact integrity** means that the bytes still match the recorded
   SHA-256 values.
2. **Environment replay** means that the same source, Python packages, base
   image, thread settings, and configuration are used.
3. **Independent reproduction** requires a separately authored implementation
   and preferably another language/toolchain.  The Python Cartesian checker is
   a deliberately separated code path, but it is not this third level.

Floating-point arrays are not promised to be byte-identical across operating
systems, NumPy wheels, BLAS implementations, CPU dispatch paths, or thread
counts.  Differences near machine precision can be numerically harmless while
still changing an NPZ checksum.  A checksum therefore authenticates one run;
it does not replace tolerance-based cross-platform validation.

## Runtime provenance v2

New artifacts use `ns-certificate-lab/runtime-provenance/v2`.  In addition to
the v1 Python, NumPy, platform, Git, and source-fingerprint fields, v2 records:

- machine, executable architecture, and byte order;
- the values (including explicit `null`) of common BLAS/OpenMP thread
  environment variables;
- a normalized JSON copy of `numpy.__config__.CONFIG`, or normalized
  `numpy.show_config()` text as a fallback;
- a SHA-256 digest of that NumPy configuration.

For normal NumPy wheels, `Build Dependencies/blas` and
`Build Dependencies/lapack` in that structure expose the linked BLAS/LAPACK
identity and build options.  Legacy runtime-provenance/v1 records remain
strictly readable, so existing candidate v2 and run-config v2 files do not need
to be rewritten.  The thread fields record environment variables, not an
independent query of an already initialized BLAS thread pool; set them before
starting Python.

The executable source fingerprint now covers:

```text
pyproject.toml
Dockerfile
.dockerignore
requirements/
src/
experiments/
configs/
scripts/
```

This is broader than the original v1-era scope.  Consequently, changing a
constraint, container definition, or release script changes the fingerprint.
Old fingerprints remain historical facts about their narrower recorded scope.

## Measured Windows reference environment

`requirements/constraints-windows-cpython311.txt` records the exact
Python-package versions measured on 2026-07-28:

- CPython 3.11.9, 64-bit;
- NumPy 1.26.4;
- pytest 9.0.3;
- pip 26.1.2 and setuptools 83.0.0;
- the exact pytest dependency versions listed in the constraints file.

Two of those pins were raised on 2026-08-02 for published security advisories:
Pygments 2.19.2 → 2.20.0 (ReDoS in the GUID regex) and setuptools 80.10.2 →
83.0.0 (MANIFEST.in exclusion bypass in sdist through a Unicode normalization
collision). Neither package enters a numerical result — Pygments is a pytest
display dependency and setuptools is build tooling — and the NumPy pin is
untouched, so the wheel observations below still describe the recorded
environment. The revised set was reinstalled from scratch on the same CPython
3.11.9 interpreter and the full suite was rerun before the change was
committed.

The installed NumPy wheel reported OpenBLAS64 0.3.23.dev, ILP64
(`USE_64BITINT=1`), OpenMP, `MAX_THREADS=2`, x86-64 little-endian, MSVC
19.29.30153, and SIMD dispatch through AVX2.  These non-Python details are
observations, not packages controlled by the constraints file; provenance v2
records them for every new run.

For a fresh PowerShell virtual environment:

```powershell
py -3.11 -m venv .venv
.\.venv\Scripts\python -m pip install "pip==26.1.2"
$env:PIP_CONSTRAINT = (Resolve-Path requirements/constraints-windows-cpython311.txt)
.\.venv\Scripts\python -m pip install -e ".[dev]"
.\.venv\Scripts\python -m pytest
Remove-Item Env:PIP_CONSTRAINT
```

These constraints are deliberately labelled for Windows CPython 3.11.  They
are not advertised as a universal lockfile and do not pin Windows itself, CPU
microcode, or the wheel repository.

## Container target

The `Dockerfile` uses CPython 3.11.9 on Debian Bookworm and pins the OCI index:

```text
python:3.11.9-slim-bookworm
sha256:8fb099199b9f2d70342674bd9dbccd3ed03a258f26bbd1d556822c6dfc60c317
```

That digest was resolved from the Docker Registry API on 2026-07-28 rather
than inferred from the tag.  The local development machine did not have a
Docker CLI, so this repository has **not yet built or executed that image**.
`requirements/constraints-container.txt` is therefore an exact intended
Python dependency target, not a verified Linux lock.  The digest prevents tag
movement, while provenance from an actual container run is still required to
record its NumPy/BLAS build and CPU dispatch.

On a machine with Docker:

```bash
docker build --tag ns-certificate-lab:repro .
docker run --rm ns-certificate-lab:repro
docker image inspect --format '{{.Id}}' ns-certificate-lab:repro
```

Record the image ID, build command, Docker version, target architecture, and
test output in release notes.  Do not claim the container is verified until
those commands pass.  The image fixes common numerical thread variables to
one; a run that overrides them will record the override in provenance.

## Review-sized release bundle

`scripts/build_release_bundle.py` creates a deterministic, stored (uncompressed)
ZIP from an allowlist.  It:

- refuses destinations outside the repository and refuses overwrite;
- rejects symlinks;
- excludes `.git`, caches, virtual environments, build directories, and
  `tests/_scratch`;
- includes source, review documents, tests, constraints, and the container
  definition;
- selects only the greatest numbered `manufactured`, `baseline`, and
  `time_convergence` output directories, plus `outputs/README.md`;
- also selects the greatest numbered `poisson_manufactured` and
  `nonlinear_control` directories when those optional evidence families are
  present;
- verifies every bundled output `.sha256` sidecar before writing;
- writes `RELEASE_INDEX.json` with every payload size and SHA-256;
- writes `SHA256SUMS` covering every payload and the index;
- uses sorted names, fixed ZIP timestamps, fixed permissions, and no
  compression so two builds from identical input bytes and Git state have
  identical archive bytes.

Build it only after source review:

```powershell
python scripts/build_release_bundle.py `
  --require-clean-source `
  --output dist/ns-singularity-certificate-lab-release.zip
```

Generated or modified files below `outputs/`, `dist/`, and `tests/_scratch/`
are ignored by the script's clean-source gate.  Other uncommitted changes are
rejected.  Independently verify the extracted archive:

```bash
cd ns-singularity-certificate-lab
sha256sum --check SHA256SUMS
```

The script does not create tags, GitHub releases, or release assets.

## Release checklist

1. Commit logically separated source and documentation changes.
2. Start from a worktree with no non-output changes and record `git rev-parse
   HEAD`.
3. Install through the applicable constraints; run the complete test suite.
4. Generate each current evidence set in a fresh output directory.  The
   experiment collects provenance before creating that directory, so it can
   truthfully record a clean pre-run source state.
5. Validate candidate/diagnostic checksums and review adverse as well as
   passing results.
6. Ensure `outputs/README.md` identifies the current generations.  If the
   generation names change, update the explicit `.gitattributes`
   `export-ignore` exceptions used by ordinary Git source archives.
7. Build the review ZIP with `--require-clean-source`; verify
   `SHA256SUMS` from a separate extraction.
8. Review `RELEASE_INDEX.json`, runtime provenance, constraints, test output,
   and any platform-dependent numerical differences.
9. Only then create an annotated tag and attach the review ZIP as a release
   asset.  Tag/release creation is an explicit maintainer action and is not
   performed by the bundle script.

Historical output directories remain in Git.  The release bundle and Git
`export-ignore` rules reduce distribution size without deleting or concealing
that evidence.

## Git line-ending migration note

The repository sets `/outputs/** -text` and `/certificates/** -text`.  This is
intentional even for JSON, CSV, SVG, and Markdown evidence: their checksum
sidecars cover exact generated bytes, so Git must not rewrite CRLF to LF during
check-in or checkout.

Some evidence predated that attribute.  Merely committing `.gitattributes`
does not rewrite an already stored Git blob.  When introducing this rule, stage
the tracked evidence again after the attribute is active:

```bash
git add --renormalize outputs certificates
```

Review those staged changes as line-ending-byte preservation, then verify
manifest and sidecar hashes from a fresh checkout or `git archive`.  A release
must not proceed if the blob bytes obtained from Git differ from the bytes
named by an included manifest.  The custom release bundle reads raw worktree
bytes and verifies sidecars plus each selected experiment `manifest.json`
files table before writing, but that does not excuse a corrupt Git source
archive.
