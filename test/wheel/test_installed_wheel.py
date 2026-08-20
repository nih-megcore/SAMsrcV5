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
from samsrcv5.fiducials import convert_json_fids_to_head
from samsrcv5.launcher import _native_path
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
    assert samsrcv5.__version__ == "5.1.0"
    assert importlib.metadata.version("samsrcv5") == "5.1.0"
    root = importlib.resources.files("samsrcv5")
    for relative in (
        "data/master+orig.HEAD",
        "data/master+orig.BRIK.gz",
        "data/ROIbuilder.ui",
        "licenses/THIRD_PARTY_NOTICES.md",
    ):
        assert root.joinpath(relative).is_file(), relative
    if os.name == "nt":
        package_dir = Path(samsrcv5.__file__).resolve().parent
        dll_dir = package_dir.parent / "samsrcv5.libs"
        assert dll_dir.is_dir()
        assert any(dll_dir.glob("*.dll"))


def test_all_console_entry_points_are_installed() -> None:
    expected = {
        "1dstats",
        "3dNormalize",
        "convert_json_fids_to_head",
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
        "sam_epi",
        "sam_ers",
        "sam_param_gui",
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
    assert callable(convert_json_fids_to_head)


def test_native_help_and_error_contracts() -> None:
    for command in ("sam_cov", "sam_wts", "sam_3d", "sam_4d", "sam_epi", "sam_ers", "sam_power"):
        result = run(command, "-h")
        assert result.returncode == 0
        assert "Usage:" in result.stdout + result.stderr
        assert "Version 5.1.0" in result.stdout + result.stderr
    result = run("sam_cov", "--CovBand", "bad", "70")
    assert result.returncode != 0
    assert "badly formed number" in result.stderr


def test_1dstats_round_trip() -> None:
    result = run("1dstats", "-q", input_text="1\n2\n3\n4\n")
    assert result.returncode == 0
    values = [float(value) for value in result.stdout.split()]
    assert values[:4] == pytest.approx([4.0, 2.5, 5.0 / 3.0, 10.0])


def test_windows_native_path_includes_repaired_dlls(tmp_path: Path) -> None:
    executable = tmp_path / "site-packages/samsrcv5/_bin/sam_wts.exe"
    executable.parent.mkdir(parents=True)
    dll_dir = tmp_path / "site-packages/samsrcv5.libs"
    dll_dir.mkdir()

    search_path = _native_path(executable, "original-path", windows=True)
    entries = search_path.split(os.pathsep)
    assert entries == [str(executable.parent), str(dll_dir), "original-path"]


def test_normalize_help_does_not_require_afni() -> None:
    result = run("3dNormalize", "--help")
    assert result.returncode == 0
    assert "Scale an AFNI or SAM volume" in result.stdout


def test_fiducial_converter_help_does_not_require_afni() -> None:
    result = run("convert_json_fids_to_head", "--help")
    assert result.returncode == 0
    assert "BIDS T1w NIfTI" in result.stdout


def test_orthohull_nifti_reports_missing_json_before_afni(tmp_path: Path) -> None:
    image = tmp_path / "subject_T1w.nii"
    image.touch()
    result = run("orthohull", str(image))
    assert result.returncode != 0
    assert "JSON sidecar does not exist" in result.stderr


def test_normalize_tail_statistics_and_compressed_view(tmp_path: Path) -> None:
    assert _tail_sd([-1.0, -3.0]) == pytest.approx((20.0 / 3.0) ** 0.5)
    prefix = tmp_path / "image+acpc"
    prefix.with_suffix(".BRIK.gz").touch()
    assert _view(str(prefix)) == ("acpc", True)
