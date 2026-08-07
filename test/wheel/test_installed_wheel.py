from __future__ import annotations

import importlib.metadata
import importlib.resources
import os
import shutil
import subprocess
import sysconfig
from pathlib import Path

import pytest
import samsrcv5
from samsrcv5.normalize import _tail_sd, _view


def run(command: str, *args: str, input_text: str | None = None):
    executable = shutil.which(command)
    if executable is None:
        candidate = Path(sysconfig.get_path("scripts")) / command
        if candidate.is_file():
            executable = str(candidate)
        elif os.name == "nt" and candidate.with_suffix(candidate.suffix + ".exe").is_file():
            executable = str(candidate.with_suffix(candidate.suffix + ".exe"))
    if executable is None and os.name == "nt" and command.endswith(".py"):
        executable = shutil.which(command.removesuffix(".py"))
    assert executable, f"console command is missing: {command}"
    return subprocess.run(
        [executable, *args],
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )


def test_metadata_and_resources() -> None:
    assert samsrcv5.__version__ == "5.0.0"
    assert importlib.metadata.version("samsrcv5") == "5.0.0"
    root = importlib.resources.files("samsrcv5")
    for relative in (
        "data/master+orig.HEAD",
        "data/master+orig.BRIK.gz",
        "data/ROIbuilder.ui",
        "licenses/THIRD_PARTY_NOTICES.md",
    ):
        assert root.joinpath(relative).is_file(), relative


def test_all_console_entry_points_are_installed() -> None:
    expected = {
        "1dstats",
        "3dNormalize",
        "FSnormals.py",
        "fiddist.py",
        "meshnorm",
        "mkGiiAtlas.py",
        "OPMsim",
        "orthohull",
        "orthohull.py",
        "plothull.py",
        "ply2fid",
        "ROIbuilder",
        "sam_3d",
        "sam_4d",
        "sam_cov",
        "sam_power",
        "sam_simulate",
        "sam_wts",
    }
    scripts = {
        entry.name
        for entry in importlib.metadata.entry_points(group="console_scripts")
        if entry.dist and entry.dist.name == "samsrcv5"
    }
    assert expected <= scripts


def test_native_help_and_error_contracts() -> None:
    for command in ("sam_cov", "sam_wts", "sam_3d", "sam_4d", "sam_power"):
        result = run(command, "-h")
        assert result.returncode == 0
        assert "Usage:" in result.stdout + result.stderr
        assert "Version 5.0" in result.stdout + result.stderr
    result = run("sam_cov", "--CovBand", "bad", "70")
    assert result.returncode != 0
    assert "badly formed number" in result.stderr


def test_1dstats_round_trip() -> None:
    result = run("1dstats", "-q", input_text="1\n2\n3\n4\n")
    assert result.returncode == 0
    values = [float(value) for value in result.stdout.split()]
    assert values[:4] == pytest.approx([4.0, 2.5, 5.0 / 3.0, 10.0])


def test_normalize_help_does_not_require_afni() -> None:
    result = run("3dNormalize", "--help")
    assert result.returncode == 0
    assert "Scale an AFNI or SAM volume" in result.stdout


def test_normalize_tail_statistics_and_compressed_view(tmp_path: Path) -> None:
    assert _tail_sd([-1.0, -3.0]) == pytest.approx((20.0 / 3.0) ** 0.5)
    prefix = tmp_path / "image+acpc"
    prefix.with_suffix(".BRIK.gz").touch()
    assert _view(str(prefix)) == ("acpc", True)
