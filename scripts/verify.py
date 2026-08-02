"""One-command verification entry point.

Runs the repository's full machine-checkable suite and then prints, in fixed
wording, exactly what that run does and does not establish.  The scope notice
is the point of this script: a green test run is easy to over-read, so the
boundary is printed whether or not the reader opens the README.

Usage::

    python scripts/verify.py            # install dependencies, then verify
    python scripts/verify.py --no-install

Exit status is the test suite's own.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

VERIFIED = """\
WHAT THIS RUN VERIFIED
  * every independent certificate checker, together with its tamper-rejection
    battery (payloads that misstate a constant, a hypothesis, or a disclaimer
    are rejected);
  * the spectral-front identities of the main theorem, evaluated in exact
    rational arithmetic on finite trigonometric fields, including the
    telescoping of the defect decomposition;
  * the coherent-family certificates that constitute the finite evidence for
    Hypothesis L*, including the exact symmetry identities and the exact value
    of K at small band sizes;
  * the exclusion certificates for the mechanisms this repository tested and
    rejected;
  * the README's own claim guards (cited theorems exist, linked paths resolve,
    the headline numbers match the artifact they cite).
"""

NOT_VERIFIED = """\
WHAT THIS RUN DID NOT VERIFY
  * Hypothesis L* itself.  It is open.  The certificates above are finite
    evidence at small band sizes plus a float continuation; they are not a
    proof, and the conditional Proposition depends on L* being true.
  * The infinite-dimensional analysis of the paper.  Those proofs are on
    paper.  They were reviewed adversarially by AI agents, which is not human
    peer review, and they are not machine-checked.
  * The Lean development.  Lean is deliberately outside this command and
    outside CI, because it needs elan/lake and a large mathlib cache.  To run
    it:  cd formal && lake exe cache get && lake build && lake env lean
    AxiomAudit.lean
  * The multi-hour Track P chain replays, which are excluded from the default
    suite for runtime reasons; docs/reproducibility.md gives their commands.

A green run above is evidence that the finite, machine-checkable layer is
internally consistent and independently recomputable.  It is not a proof of
any statement about the Navier-Stokes equations, and this repository does not
claim a solution to the Clay problem in either direction.
"""


def _run(command: list[str], *, why: str) -> int:
    print(f"\n$ {' '.join(command)}\n  ({why})", flush=True)
    return subprocess.run(command, cwd=ROOT).returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-install",
        action="store_true",
        help="skip the editable install and use the current environment",
    )
    arguments = parser.parse_args()

    if sys.version_info < (3, 10):
        print(f"Python 3.10 or newer is required; this is {sys.version.split()[0]}")
        return 1

    if not arguments.no_install:
        code = _run(
            [sys.executable, "-m", "pip", "install", "-e", ".[dev]", "-q"],
            why="numpy and pytest only; scipy is optional and one test skips without it",
        )
        if code != 0:
            print("\nInstallation failed; rerun with --no-install to use this "
                  "environment as it stands.")
            return code

    code = _run([sys.executable, "-m", "pytest", "-q"], why="the full suite")

    print()
    print("=" * 72)
    print(VERIFIED)
    print(NOT_VERIFIED, end="")
    print("=" * 72)
    if code != 0:
        print("\nThe suite did NOT pass; the scope notice above describes what a "
              "passing run would have covered.")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
