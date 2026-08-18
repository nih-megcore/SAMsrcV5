"""Console entry points for bundled SAM executables and legacy Python tools."""

from __future__ import annotations

import os
import runpy
import subprocess
import sys
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
    _legacy("orthohull.py")


def orthohull_python() -> None:
    _legacy("orthohull.py")


def plothull() -> None:
    _legacy("plothull.py")


def ply2fid() -> None:
    _legacy("ply2fid.py")


def roi_builder() -> None:
    _legacy("ROIbuilder.py")
