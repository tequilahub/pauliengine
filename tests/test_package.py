from __future__ import annotations

import importlib.metadata

import pauliengine as m


def test_version():
    assert importlib.metadata.version("pauliengine") == m.__version__
