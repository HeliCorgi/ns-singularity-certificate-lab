"""Semantic staleness guards for README.md and the verified-results registry.

These tests are deliberately *semantic* rather than wording-pinned: they check
that what the README asserts still matches the tree, not that particular
sentences survive editing.  What they enforce:

- every Lean theorem the README cites exists in ``formal/``;
- every repository path the README links resolves;
- the headline representative result's numbers match the artifact it cites,
  and the pinned reproduction commit really contains that artifact;
- the promotion/hybrid block appears exactly when every audited external
  record carries ``proved: True`` with complete closure metadata;
- nothing whose ``lean_formalised`` flag is false is described as Lean-only,
  and nothing conditional is presented as unconditional;
- deleted stale claims do not reappear;
- the branch story is consistent with the branches that exist;
- quick-verification commands point at files that exist;
- no LaTeX is broken and no GitHub-incompatible math delimiter is used.
"""

from __future__ import annotations

import json
import re
import subprocess
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
README = (ROOT / "README.md").read_text(encoding="utf-8")
REGISTRY_PATH = ROOT / "docs" / "verified_results.md"

SECTION_RE = re.compile(r"<!-- MCR:BEGIN.*?<!-- MCR:END -->", re.DOTALL)
PROMOTION_RE = re.compile(r"<!-- PROMOTION:BEGIN.*?<!-- PROMOTION:END -->", re.DOTALL)

#: The representative result the README leads with.
REPRESENTATIVE_ARTIFACT = (
    "outputs/track_p_chain_reissued_v2/"
    "reissued_h3chain_strict_same_step_P1_nu_1over10.json"
)


@pytest.fixture(scope="module")
def section() -> str:
    match = SECTION_RE.search(README)
    assert match, "README must contain the MCR:BEGIN/END machine-checkable section"
    return match.group(0)


@pytest.fixture(scope="module")
def lean_sources() -> str:
    parts = [
        path.read_text(encoding="utf-8")
        for path in sorted((ROOT / "formal").rglob("*.lean"))
        if ".lake" not in path.parts
    ]
    assert parts, "formal/ must contain Lean sources"
    return "\n".join(parts)


@pytest.fixture(scope="module")
def flat() -> str:
    """README with soft line wrapping collapsed, for phrase checks."""
    return re.sub(r"\s+", " ", README)


@pytest.fixture(scope="module")
def asserted() -> str:
    """README minus the block that *disclaims* phrases rather than using them.

    The vocabulary section lists the overclaims this project refuses to make.
    A naive substring scan would fire on the disclaimer itself, so the block
    is delimited and removed before scanning for asserted overclaims.
    """
    stripped = re.sub(
        r"<!-- NOTCLAIMED:BEGIN.*?<!-- NOTCLAIMED:END -->",
        "",
        README,
        flags=re.DOTALL,
    )
    return re.sub(r"\s+", " ", stripped)


# --------------------------------------------------------------------------- #
# citations resolve                                                            #
# --------------------------------------------------------------------------- #


def test_cited_identifiers_exist(lean_sources):
    """Cited Lean theorems and Python checker entry points must both exist.

    A dotted name whose prefix is a module of ``ns_certificate_lab`` is checked
    as a Python attribute; everything else is checked against ``formal/``.
    """
    cited = set(re.findall(r"`([A-Za-z_][A-Za-z0-9_.']*)`", README))
    ignore = {
        "sorry", "admit", "axiom", "main", "fable5-mainline",
        "lean-toolchain", "lake-manifest.json", "C_kato", "C_shift",
        "paper_proof_status", "lean_formalised", "external_theorem_status",
        "audited_not_formalised", "slab_budget_exhausted", "proved",
        "snapshot_certificate", "control_linear_coefficient",
    }
    src = ROOT / "src" / "ns_certificate_lab"
    missing_lean, missing_python = [], []
    for name in cited:
        if name in ignore or name.endswith((".py", ".json", ".md", ".lean", ".toml")):
            continue
        if "_" not in name and "." not in name:
            continue
        module, _, attribute = name.rpartition(".")
        if module and (src / f"{module}.py").exists():
            text = (src / f"{module}.py").read_text(encoding="utf-8")
            defined = re.search(
                rf"^(def |class ){re.escape(attribute)}\b"
                rf"|^{re.escape(attribute)}\s*[:=]",
                text,
                re.MULTILINE,
            )
            if not defined:
                missing_python.append(name)
            continue
        base = name.split(".")[-1]
        if not re.search(rf"\b{re.escape(base)}\b", lean_sources):
            missing_lean.append(name)
    assert not missing_lean, f"README cites Lean names absent from formal/: {missing_lean}"
    assert not missing_python, (
        f"README cites checker entry points that do not exist: {missing_python}"
    )


def test_linked_paths_exist():
    links = re.findall(r"\]\(([^)#\s]+)\)", README)
    missing = [
        link
        for link in links
        if not link.startswith("http") and not (ROOT / link).exists()
    ]
    assert not missing, f"README links to nonexistent paths: {missing}"


def test_quick_verification_targets_exist():
    """The files the quick-verification block drives must be present."""
    for path in ("pyproject.toml", "tests/test_readme_claims.py",
                 "formal/AxiomAudit.lean", "formal/lean-toolchain",
                 "formal/lake-manifest.json"):
        assert (ROOT / path).exists(), f"quick verification references {path}"
    assert "python -m pytest" in README
    assert "lake build" in README
    assert "lake env lean AxiomAudit.lean" in README


# --------------------------------------------------------------------------- #
# the headline result matches its artifact                                     #
# --------------------------------------------------------------------------- #


def test_representative_result_matches_its_artifact(flat):
    artifact = ROOT / REPRESENTATIVE_ARTIFACT
    assert artifact.exists(), f"the representative artifact {artifact} is missing"
    payload = json.loads(artifact.read_text(encoding="utf-8"))
    order = payload["norm_order"]

    assert REPRESENTATIVE_ARTIFACT in README, (
        "the README must link the artifact backing its representative result"
    )
    # certified interval, exactly as stored
    assert payload["certified_final_time"] in flat, (
        f"README does not state the artifact's certified interval "
        f"{payload['certified_final_time']}"
    )
    # tube radius, to the precision the README quotes
    radius = float(Fraction(payload[f"max_radius_h{order}_dot"]))
    assert f"{radius:.6f}" in flat, (
        f"README does not state the artifact's tube radius {radius:.6f}"
    )
    # Sobolev order and slab count
    assert f"H^{order}" in flat or f"H^{{{order}}}" in flat or f"\\dot H^{order}" in flat
    assert str(payload["slab_count"]) in flat
    # viscosity
    assert payload["viscosity"] in flat
    # the stop reason must be reported, not hidden
    assert payload["termination"]["reason"] in README


def test_pinned_commit_exists_and_contains_the_artifact(flat):
    """The commit the README pins for exact reproduction must be real."""
    match = re.search(r"git checkout ([0-9a-f]{7,40})", README)
    assert match, "README must pin a commit for exact reproduction"
    commit = match.group(1)
    listing = subprocess.run(
        ["git", "cat-file", "-e", f"{commit}^{{commit}}"],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert listing.returncode == 0, f"pinned commit {commit} does not exist"
    tree = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", commit, REPRESENTATIVE_ARTIFACT],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert REPRESENTATIVE_ARTIFACT in tree.stdout, (
        f"pinned commit {commit} does not contain {REPRESENTATIVE_ARTIFACT}"
    )


# --------------------------------------------------------------------------- #
# classification integrity                                                     #
# --------------------------------------------------------------------------- #


def test_promotion_block_matches_recorded_ext_state(section):
    from ns_certificate_lab.torus_aposteriori import EXTERNAL_THEOREMS_AUDITED

    all_closed = all(
        block.get("proved") is True
        and block.get("closure", {}).get("method") == "audited_paper_proof"
        and block.get("closure", {}).get("lean_formalised") is False
        and block.get("closure", {}).get("axiomatised_in_lean") is False
        for block in EXTERNAL_THEOREMS_AUDITED.values()
    )
    promotion = PROMOTION_RE.search(section)
    if promotion:
        assert all_closed, (
            "README carries the hybrid/promotion block while an audited EXT "
            "record is not fully closed"
        )
        text = re.sub(r"\s+", " ", promotion.group(0))
        assert "Lean 形式化済みを意味しません" in text, (
            "the hybrid block must state that proved:true is not a Lean "
            "formalisation"
        )
        assert "監査済み" in text


def test_audited_but_unformalised_is_never_called_lean_only(asserted):
    """Nothing with lean_formalised false may be presented as a Lean result."""
    from ns_certificate_lab.torus_aposteriori import EXTERNAL_THEOREMS_AUDITED

    unformalised = [
        name
        for name, block in EXTERNAL_THEOREMS_AUDITED.items()
        if block.get("closure", {}).get("lean_formalised") is False
    ]
    assert unformalised, "the audited records should mark the analysis unformalised"
    forbidden = (
        "Lean だけで周期 Navier–Stokes",
        "fully Lean-verified PDE theorem",
        "Lean-verified PDE theorem",
        "Lean で証明された周期 Navier–Stokes 解の存在",
    )
    for phrase in forbidden:
        assert phrase not in asserted, f"README claims too much: {phrase!r}"


def test_conditional_artifacts_are_not_presented_as_unconditional(section, flat):
    """v1 payloads keep the conditional statement; the README must say so."""
    summary_path = ROOT / "outputs" / "track_p_chain_v1" / "summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        assert "conditional" in summary.get("disclaimer", "")
    assert "C1." in section, "the conditional-certificate section must exist"
    assert "EXT 仮定に" in flat and "条件付き" in flat, (
        "section C1 must state that the v1 PDE conclusions stay conditional"
    )


def test_vocabulary_definitions_present(section):
    for term in ("Lean-verified", "certificate-verified", "hybrid",
                 "conditional", "numerically observed"):
        assert term in section, f"vocabulary entry missing: {term}"


def test_forbidden_overclaims_absent(asserted):
    for phrase in ("Clay 問題への部分解", "大域正則性証明", "特異点証明を与える",
                   "fully verified", "machine-checked PDE theorem"):
        assert phrase not in asserted, f"forbidden claim present: {phrase!r}"


def test_open_dependencies_are_listed(section):
    for item in ("G-DINI", "HS-5", "NT-N1", "H3"):
        assert item in section, f"open dependency not listed: {item}"


def test_auditor_description_is_accurate(flat):
    """The audits were agent passes; the README must not imply human review."""
    if "監査" in flat:
        assert "エージェント監査" in flat or "independent agent review" in flat, (
            "the README must say the audit passes were agent-run, not human "
            "peer review"
        )
        assert "人間の査読ではありません" in flat or "not human peer review" in flat


# --------------------------------------------------------------------------- #
# stale content must not come back                                             #
# --------------------------------------------------------------------------- #


def test_deleted_stale_claims_do_not_reappear(flat):
    stale = (
        "既知に滑らかな人工場と",          # pre-Track-P scope sentence
        "非特異な減衰対照だけです",
        "区間演算、厳密な打切り誤差評価",   # "not yet implemented" list
        "圧力回復の独立実装はまだありません",
        "候補探索は未着手",
        "全空間solverは未実装",
    )
    for phrase in stale:
        assert phrase not in flat, f"stale claim reappeared: {phrase!r}"


def test_branch_story_is_consistent():
    """Branches the README names must exist, with main as the stable one."""
    result = subprocess.run(
        ["git", "branch", "-a", "--format=%(refname:short)"],
        cwd=ROOT, capture_output=True, text=True,
    )
    refs = set(result.stdout.split())
    for branch in ("main", "fable5-mainline"):
        assert branch in README, f"README must describe the {branch} branch"
        assert any(branch == r or r.endswith(f"/{branch}") for r in refs), (
            f"README names branch {branch}, which does not exist"
        )
    assert "最新の安定版" in README and "開発版" in README


# --------------------------------------------------------------------------- #
# rendering                                                                    #
# --------------------------------------------------------------------------- #


def test_no_broken_latex():
    """Catch backslash-eaten commands like 'abla' or 'angle' from '\\nabla'."""
    for broken in (r"\cdot abla", " abla)", "angle=0", r"\r angle",
                   "mathbb R^3\\)"):
        assert broken not in README, f"broken LaTeX in README: {broken!r}"
    # a bare 'abla' or 'rangle' at line start is the classic corruption
    for line in README.splitlines():
        stripped = line.strip()
        assert not stripped.startswith(("abla", "angle", "rangle", "nabla ")), (
            f"line looks like corrupted LaTeX: {line!r}"
        )


def test_math_uses_github_compatible_delimiters():
    r"""GitHub renders $...$ and $$...$$, never \( \) or \[ \]."""
    assert r"\(" not in README, (
        r"README uses \( \) math delimiters, which GitHub does not render"
    )
    assert r"\[" not in README, (
        r"README uses \[ \] math delimiters, which GitHub does not render"
    )


def test_no_duplicate_section_numbers():
    headings = re.findall(r"^#{2,3}\s+(\d+(?:\.\d+)?)\.", README, re.MULTILINE)
    duplicates = {h for h in headings if headings.count(h) > 1}
    assert not duplicates, f"duplicate section numbers in README: {duplicates}"


def test_no_hardcoded_drifting_counts(flat):
    assert "変動" in flat, "the README must note that counts drift"
    assert not re.search(r"\b1[0-9]{3} passed\b", README), (
        "do not hard-code pytest counts"
    )
    assert not re.search(r"全\s*\d{3}\s*定理", README), (
        "do not hard-code the audited theorem count; point at AxiomAudit.lean"
    )


# --------------------------------------------------------------------------- #
# registry agreement                                                           #
# --------------------------------------------------------------------------- #


def test_registry_exists_and_is_linked():
    assert REGISTRY_PATH.exists(), "docs/verified_results.md must exist"
    assert "docs/verified_results.md" in README, (
        "the README must point at the registry for the complete list"
    )


def test_registry_paths_resolve():
    registry = REGISTRY_PATH.read_text(encoding="utf-8")
    links = re.findall(r"\]\(([^)#\s]+)\)", registry)
    missing = []
    for link in links:
        if link.startswith("http"):
            continue
        if not (REGISTRY_PATH.parent / link).resolve().exists() and not (
            ROOT / link
        ).exists():
            missing.append(link)
    assert not missing, f"registry links to nonexistent paths: {missing}"


def test_registry_covers_the_representative_result():
    registry = REGISTRY_PATH.read_text(encoding="utf-8")
    assert "track_p_chain_reissued_v2" in registry, (
        "the registry must carry the reissued chain family the README "
        "leads with"
    )


def test_no_sorry_admit_axiom_in_lean_sources(lean_sources):
    code = re.sub(r"/-.*?-/", "", lean_sources, flags=re.DOTALL)
    code = re.sub(r"--[^\n]*", "", code)
    for needle in ("sorry", "admit"):
        assert not re.search(rf"(^|\s){needle}\b", code, re.MULTILINE), (
            f"Lean sources contain {needle}"
        )
    assert not re.search(r"^\s*axiom\s", code, re.MULTILINE), (
        "Lean sources declare a project axiom"
    )


def test_ci_coverage_claim_matches_the_workflow(flat):
    """The README says CI re-runs Python only; verify that stays true.

    If a Lean build/audit job is ever added to CI, this fires so the README
    stops under-claiming (and vice versa if the Python job is removed).
    """
    workflows = sorted((ROOT / ".github" / "workflows").glob("*.yml"))
    assert workflows, "CI workflow directory must exist"
    text = "\n".join(p.read_text(encoding="utf-8") for p in workflows)
    runs_python = "pytest" in text
    runs_lean = bool(re.search(r"^\s*(-\s*)?run:.*\b(lake|lean)\b", text, re.M))
    assert runs_python, "CI no longer runs pytest; the README says it does"
    if runs_lean:
        pytest.fail(
            "CI now runs Lean, but the README states the Lean build and axiom "
            "audit are not part of CI — update the README"
        )
    assert "Lean の build と公理監査は CI に含まれていません" in flat, (
        "the README must state honestly that Lean is not covered by CI"
    )


def test_checker_scope_is_not_overclaimed(section):
    """Section B must not promise recomputation the slab checker disclaims."""
    from ns_certificate_lab import torus_aposteriori

    doc = torus_aposteriori.verify_torus_slab_certificate.__doc__ or ""
    if "not" in doc and "re-run the interval convolutions" in doc:
        assert "checker は再計算しない" in section or "builder 側" in section, (
            "the slab checker explicitly does NOT re-run the M_j interval "
            "convolutions; section B must say so"
        )


def test_every_vocabulary_term_has_a_representative_section(section):
    """All five statuses must be represented, not just the flattering four."""
    for heading in ("### A.", "### B.", "### C1.", "### C2.", "### C3.",
                    "### C4.", "### D."):
        assert heading in section, f"missing section {heading}"
    assert "numerically observed" in section.lower() or "浮動小数点による観測" in section


def test_enforced_contract_strings_are_quoted_verbatim(flat):
    """C1/C2 must quote the checker-enforced wording, not paraphrase it."""
    from ns_certificate_lab.torus_chain import ALLOWED_WORDING, AUDITED_KIND

    unquoted = re.sub(r"^>\s?", "", README, flags=re.MULTILINE)
    flat_readme = re.sub(r"\s+", " ", unquoted)
    for sentence in ALLOWED_WORDING:
        assert re.sub(r"\s+", " ", sentence) in flat_readme, (
            f"README must quote the enforced wording verbatim: {sentence[:60]!r}"
        )
    head = re.sub(r"\s+", " ", AUDITED_KIND)[:60]
    assert head in flat_readme, "README must quote AUDITED_KIND verbatim"
