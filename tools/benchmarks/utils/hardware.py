"""Collect hardware and software fingerprint for a benchmark run.

The goal is that a saved result JSON is self-describing: a plot generated
from it can be labeled with the machine it was produced on, and two result
files can be compared safely (or flagged as incomparable).
"""

from __future__ import annotations

import platform
import subprocess
import sys
from typing import Any


def _try(func):
    try:
        return func()
    except Exception:  # noqa: BLE001 - hardware probing is best-effort
        return None


def _cpu_model() -> str | None:
    system = platform.system()
    if system == "Windows":
        out = _try(lambda: subprocess.check_output(
            ["wmic", "cpu", "get", "Name"], text=True, stderr=subprocess.DEVNULL,
        ))
        if out:
            lines = [line.strip() for line in out.splitlines() if line.strip()]
            if len(lines) >= 2:
                return lines[1]
        return platform.processor() or None
    if system == "Darwin":
        return _try(lambda: subprocess.check_output(
            ["sysctl", "-n", "machdep.cpu.brand_string"], text=True,
        ).strip())
    if system == "Linux":
        def read_cpuinfo() -> str | None:
            with open("/proc/cpuinfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        return line.split(":", 1)[1].strip()
            return None
        return _try(read_cpuinfo) or platform.processor() or None
    return platform.processor() or None


def _total_ram_bytes() -> int | None:
    system = platform.system()
    if system == "Linux":
        def read_meminfo() -> int:
            with open("/proc/meminfo", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        return kb * 1024
            raise RuntimeError("MemTotal not found")
        return _try(read_meminfo)
    if system == "Darwin":
        out = _try(lambda: subprocess.check_output(["sysctl", "-n", "hw.memsize"], text=True))
        return int(out.strip()) if out else None
    if system == "Windows":
        out = _try(lambda: subprocess.check_output(
            ["wmic", "ComputerSystem", "get", "TotalPhysicalMemory"],
            text=True, stderr=subprocess.DEVNULL,
        ))
        if out:
            lines = [line.strip() for line in out.splitlines() if line.strip()]
            if len(lines) >= 2 and lines[1].isdigit():
                return int(lines[1])
    return None


def _cpu_count_logical() -> int | None:
    import os
    return _try(os.cpu_count)


def _pkg_version(name: str) -> str | None:
    try:
        from importlib.metadata import version
        return version(name)
    except Exception:  # noqa: BLE001
        return None


def _git_commit() -> str | None:
    out = _try(lambda: subprocess.check_output(
        ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL,
    ))
    return out.strip() if out else None


def collect() -> dict[str, Any]:
    """Return a JSON-serializable dict describing the host and package versions."""
    return {
        "os": {
            "system": platform.system(),
            "release": platform.release(),
            "version": platform.version(),
            "machine": platform.machine(),
        },
        "cpu": {
            "model": _cpu_model(),
            "logical_cores": _cpu_count_logical(),
        },
        "memory": {
            "total_bytes": _total_ram_bytes(),
        },
        "python": {
            "version": sys.version.split()[0],
            "implementation": platform.python_implementation(),
        },
        "packages": {
            "pauliengine": _pkg_version("pauliengine"),
            "numpy": _pkg_version("numpy"),
        },
        "git": {
            "commit": _git_commit(),
        },
    }


def short_label(hw: dict[str, Any]) -> str:
    """One-line label to stamp on plots so results aren't mixed up."""
    cpu = hw.get("cpu", {}).get("model") or "unknown CPU"
    os_name = hw.get("os", {}).get("system") or "?"
    ver = hw.get("packages", {}).get("pauliengine") or "?"
    return f"{cpu} | {os_name} | pauliengine={ver}"
