"""Test fixtures that keep every generated file inside the repository."""

from __future__ import annotations

from pathlib import Path
import shutil
import uuid

import pytest


@pytest.fixture
def scratch_dir() -> Path:
    root = (Path(__file__).resolve().parent / "_scratch").resolve()
    destination = (root / uuid.uuid4().hex).resolve()
    if root not in destination.parents:
        raise RuntimeError("refusing to create a test directory outside tests/_scratch")
    destination.mkdir(parents=True)
    try:
        yield destination
    finally:
        if destination.exists():
            resolved = destination.resolve()
            if root not in resolved.parents:
                raise RuntimeError("refusing to remove a directory outside tests/_scratch")
            shutil.rmtree(resolved)
        if root.exists() and not any(root.iterdir()):
            root.rmdir()
