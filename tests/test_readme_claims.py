"""Staleness guards for the README machine-checkable section and registry.

Runs in every CI push.  What it enforces:

- every Lean theorem name cited in the README section and the registry
  exists in the ``formal/`` sources;
- every repository path (files/dirs) linked in the section exists;
- the EXT promotion block appears ONLY while every record of
  ``EXTERNAL_THEOREMS_AUDITED`` carries ``proved: True`` with complete
  closure metadata (and disappears if any flag regresses);
- the promotion block carries the mandatory scope caveats;
- v1 conditional artifacts are not presented as unconditional (their
  summaries keep the conditional disclaimer);
- no ``sorry``/``admit``/project ``axiom`` appears in the Lean sources;
- drifting counts are not hard-coded as fixed facts in the section.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
README = (ROOT / "README.md").read_text(encoding="utf-8")
REGISTRY_PATH = ROOT / "docs" / "verified_results.md"

SECTION_RE = re.compile(r"<!-- MCR:BEGIN.*?<!-- MCR:END -->", re.DOTALL)
PROMOTION_RE = re.compile(r"<!-- PROMOTION:BEGIN.*?<!-- PROMOTION:END -->", re.DOTALL)


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


def test_cited_lean_theorems_exist(section, lean_sources):
    cited = set(re.findall(r"`([A-Za-z_][A-Za-z0-9_.']*)`", section))
    lean_like = {
        name
        for name in cited
        if re.match(r"^[a-zA-Z]", name)
        and ("_" in name or "." in name)
        and not name.endswith((".py", ".json", ".md", ".lean"))
        and name
        not in {
            "sorry", "admit", "axiom", "lake env lean AxiomAudit.lean",
            "proved: true", "lean_formalised: false", "C_kato", "C_shift",
            "lean-toolchain", "lake-manifest.json", "sorry`/`admit`/`axiom",
        }
    }
    assert lean_like, "the README section should cite Lean theorem names"
    missing = []
    for name in lean_like:
        base = name.split(".")[-1] if "." in name else name
        if not re.search(rf"\b{re.escape(base)}\b", lean_sources):
            missing.append(name)
    assert not missing, f"README cites Lean names not found in formal/: {missing}"


def test_linked_paths_exist(section):
    links = re.findall(r"\]\(([^)#\s]+)\)", section)
    missing = [
        link
        for link in links
        if not link.startswith("http") and not (ROOT / link).exists()
    ]
    assert not missing, f"README links to nonexistent paths: {missing}"


def test_registry_exists_and_paths_resolve():
    assert REGISTRY_PATH.exists(), "docs/verified_results.md must exist"
    registry = REGISTRY_PATH.read_text(encoding="utf-8")
    links = re.findall(r"\]\((\.\./[^)#\s]+|[^)#\s:]+)\)", registry)
    missing = []
    for link in links:
        if link.startswith("http"):
            continue
        candidate = (REGISTRY_PATH.parent / link).resolve()
        alternative = (ROOT / link).resolve()
        if not candidate.exists() and not alternative.exists():
            missing.append(link)
    assert not missing, f"registry links to nonexistent paths: {missing}"


def test_promotion_block_consistent_with_ext_flags(section):
    from ns_certificate_lab.torus_aposteriori import EXTERNAL_THEOREMS_AUDITED

    promotion = PROMOTION_RE.search(section)
    all_closed = all(
        block.get("proved") is True
        and block.get("closure", {}).get("method") == "audited_paper_proof"
        and block.get("closure", {}).get("lean_formalised") is False
        and block.get("closure", {}).get("axiomatised_in_lean") is False
        for block in EXTERNAL_THEOREMS_AUDITED.values()
    )
    if promotion:
        assert all_closed, (
            "README carries the promotion block while an audited EXT record "
            "is not fully closed — remove the block or restore the closure"
        )
        text = re.sub(r"\s+", " ", promotion.group(0))
        for caveat in ("大域正則性", "特異点構成ではな", "Clay 問題の解決ではない",
                       "証明手法の到達範囲"):
            assert caveat in text, f"promotion block must carry caveat: {caveat}"
        assert "監査済み" in text and (
            "Lean 形式化ではない" in text
            or "modulo the audited classical theorems" in text
        ), "promotion block must state the audited-paper-proof trust basis"
    else:
        # With everything closed the promotion block is optional; with
        # anything open it must be absent — nothing further to check here.
        pass


def test_readme_never_claims_unconditional_pde_without_qualifier(section):
    from ns_certificate_lab.torus_aposteriori import EXTERNAL_THEOREMS_AUDITED

    any_open = any(
        block.get("proved") is not True
        for block in EXTERNAL_THEOREMS_AUDITED.values()
    )
    if any_open:
        assert "無条件" not in re.sub(
            r"無条件に検証できる結果|無条件外部仮定", "", section.split("### C.")[1]
        ), (
            "an EXT record is open but README section C presents the PDE "
            "tube as unconditional"
        )


def test_v1_artifacts_stay_conditional():
    summary_path = ROOT / "outputs" / "track_p_chain_v1" / "summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        assert "conditional" in summary.get("disclaimer", ""), (
            "the v1 chain summary must keep its conditional disclaimer"
        )


def test_conditional_wording_present(section):
    required = (
        "有限次元軌道、残差、control ODE、スラブ\n  連結は機械検査済み"
    )
    flat = section.replace("\n  ", " ").replace("\n", " ")
    assert "有限次元軌道、残差、control ODE、スラブ 連結は機械検査済み" in flat or \
        "有限次元軌道、残差、control ODE、スラブ連結は機械検査済み" in flat, (
        "section C must carry the mandated conditional sentence for v1"
    )
    assert "EXT-P1/P2/P3 に条件付き" in flat


def test_no_sorry_admit_axiom_in_lean_sources(lean_sources):
    # Strip comments first: doc/comment blocks legitimately DISCUSS the words
    # ("axiom is forbidden...") without declaring anything.
    code = re.sub(r"/-.*?-/", "", lean_sources, flags=re.DOTALL)
    code = re.sub(r"--[^\n]*", "", code)
    for needle in ("sorry", "admit"):
        pattern = re.compile(rf"(^|\s){needle}\b", re.MULTILINE)
        assert not pattern.search(code), f"Lean sources contain {needle}"
    axiom_pattern = re.compile(r"^\s*axiom\s", re.MULTILINE)
    assert not axiom_pattern.search(code), (
        "Lean sources declare a project axiom"
    )


def test_no_hardcoded_drifting_counts(section):
    quick = section.split("Quick verification")[1].split("Full verification")[0]
    assert "変動します" in section or "変動し" in section, (
        "the section must state that test/job counts drift"
    )
    assert not re.search(r"\b1[0-9]{3} passed\b", quick), (
        "do not hard-code test counts in the reproduce block"
    )


def test_terminology_rules_present(section):
    for term in ("Lean-verified", "certificate-verified",
                 "conditional PDE certificate", "numerically observed",
                 "candidate"):
        assert term in section, f"terminology rule missing: {term}"
