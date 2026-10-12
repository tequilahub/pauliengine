"""Sphinx configuration for the PauliEngine documentation.

Build (from the repo root, with pauliengine installed and the docs group):

    pip install . --group docs        # or: pip install sphinx breathe furo
    sphinx-build -b html docs docs/_build/html

The Python API is read from the installed ``pauliengine`` package (docstrings of
the nanobind bindings and the Python wrappers). The C++ API is read from the
headers via Doxygen XML, rendered by Breathe; Doxygen runs automatically when it
is on PATH, otherwise the C++ pages are skipped with a warning.
"""

from __future__ import annotations

import shutil
import subprocess
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

DOCS_DIR = Path(__file__).resolve().parent
DOXYGEN_XML = DOCS_DIR / "_build" / "doxygen" / "xml"

# -- Project -------------------------------------------------------------------

project = "PauliEngine"
author = "PauliEngine developers"
copyright = "PauliEngine developers"
try:
    release = version("pauliengine")
except PackageNotFoundError:
    release = "unknown"
version = ".".join(release.split(".")[:2])

# -- General -------------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "breathe",
]

exclude_patterns = ["_build"]
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "special-members": "__init__, __mul__, __add__, __sub__, __eq__, __hash__, __contains__, __len__",
    "show-inheritance": False,
}
autodoc_member_order = "groupwise"
autodoc_typehints = "description"

# -- C++ API (Doxygen XML via Breathe) -------------------------------------------


def _run_doxygen() -> bool:
    if shutil.which("doxygen") is None:
        print("warning: doxygen not found on PATH, skipping the C++ API reference")
        return False
    (DOCS_DIR / "_build" / "doxygen").mkdir(parents=True, exist_ok=True)
    subprocess.run(["doxygen", "Doxyfile"], cwd=DOCS_DIR, check=True)
    return True


HAVE_CPP_API = _run_doxygen()
breathe_projects = {"pauliengine": str(DOXYGEN_XML)}
breathe_default_project = "pauliengine"
breathe_default_members = ("members", "undoc-members")
if not HAVE_CPP_API:
    exclude_patterns.append("cpp_api.rst")

# -- HTML ------------------------------------------------------------------------

try:
    import furo  # noqa: F401

    html_theme = "furo"
except ImportError:
    html_theme = "alabaster"
html_title = f"PauliEngine {release}"
