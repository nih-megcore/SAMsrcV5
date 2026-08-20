"""Console entry points for bundled SAM executables and legacy Python tools."""

from __future__ import annotations

import getopt
import os
import runpy
import subprocess
import sys
import tempfile
from importlib.resources import as_file, files
from pathlib import Path
from typing import NoReturn


def _resource_path(*parts: str):
    return files("samsrcv5").joinpath(*parts)


def _native_path(executable: Path, current_path: str, *, windows: bool) -> str:
    native_paths = [str(executable.parent)]
    if windows:
        # delvewheel places the MinGW/OpenMP runtime DLLs beside the package
        # directory. Its import patch affects this Python process, but DLL
        # lookup for the executable child relies on PATH.
        dll_dir = executable.parent.parent.parent / "samsrcv5.libs"
        if dll_dir.is_dir():
            native_paths.append(str(dll_dir))
    if current_path:
        native_paths.append(current_path)
    return os.pathsep.join(native_paths)


def _native(name: str) -> NoReturn:
    filename = f"{name}.exe" if os.name == "nt" else name
    resource = _resource_path("_bin", filename)
    with as_file(resource) as executable:
        executable = Path(executable)
        if not executable.is_file():
            raise SystemExit(f"samsrcv5 installation is missing {filename}")
        argv = [str(executable), *sys.argv[1:]]
        env = os.environ.copy()
        env["PATH"] = _native_path(
            executable, env.get("PATH", ""), windows=os.name == "nt"
        )
        env.setdefault("HOME", str(Path.home()))
        if os.name != "nt":
            os.execve(executable, argv, env)
        raise SystemExit(subprocess.call(argv, env=env))


def _legacy(name: str) -> None:
    script = _resource_path("_legacy", name)
    with as_file(script) as script_path:
        legacy_dir = script_path.parent
        sys.path.insert(0, str(legacy_dir))
        old_argv = sys.argv
        sys.argv = [old_argv[0], *old_argv[1:]]
        try:
            runpy.run_path(str(script_path), run_name="__main__")
        finally:
            sys.argv = old_argv
            try:
                sys.path.remove(str(legacy_dir))
            except ValueError:
                pass


def _orthohull_nifti_argument(argv: list[str]) -> tuple[int, Path] | None:
    try:
        _options, arguments = getopt.getopt(argv[1:], "cqmtop:i:x:")
    except getopt.GetoptError:
        return None
    if len(arguments) not in {1, 2}:
        return None
    image = Path(arguments[0])
    if not (image.name.endswith(".nii") or image.name.endswith(".nii.gz")):
        return None
    return len(argv) - len(arguments), image


def _run_orthohull() -> None:
    nifti_argument = _orthohull_nifti_argument(sys.argv)
    if nifti_argument is None:
        _legacy("orthohull.py")
        return

    from .fiducials import FiducialConversionError, convert_json_fids_to_head

    argument_index, image = nifti_argument
    original_argv = sys.argv
    with tempfile.TemporaryDirectory(prefix="samsrcv5-orthohull-") as temporary:
        try:
            _brik, head = convert_json_fids_to_head(image, temporary)
        except FiducialConversionError as error:
            raise SystemExit(f"orthohull: error: {error}") from error

        converted_argv = original_argv.copy()
        converted_argv[argument_index] = str(head)
        sys.argv = converted_argv
        try:
            _legacy("orthohull.py")
        finally:
            sys.argv = original_argv


def stats_1d() -> NoReturn:
    _native("1dstats")


def meshnorm() -> NoReturn:
    _native("meshnorm")


def opm_sim() -> NoReturn:
    _native("OPMsim")


def sam_3d() -> NoReturn:
    _native("sam_3d")


def sam_4d() -> NoReturn:
    _native("sam_4d")


def sam_cov() -> NoReturn:
    _native("sam_cov")


def sam_ers() -> NoReturn:
    _native("sam_ers")


def sam_epi() -> NoReturn:
    _native("sam_epi")


def sam_power() -> NoReturn:
    _native("sam_power")


def sam_simulate() -> NoReturn:
    _native("sam_simulate")


def sam_wts() -> NoReturn:
    _native("sam_wts")


def fsnormals() -> None:
    _legacy("FSnormals.py")


def fiddist() -> None:
    _legacy("fiddist.py")


def mk_gii_atlas() -> None:
    _legacy("mkGiiAtlas.py")


def orthohull() -> None:
    _run_orthohull()


def orthohull_python() -> None:
    _run_orthohull()


def plothull() -> None:
    _legacy("plothull.py")


def ply2fid() -> None:
    _legacy("ply2fid.py")


def roi_builder() -> None:
    _legacy("ROIbuilder.py")
