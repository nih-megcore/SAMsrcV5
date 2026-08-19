from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import nibabel as nib
import numpy as np
import pytest

from samsrcv5 import fiducials


pytestmark = pytest.mark.unit


def write_bids_t1w(tmp_path: Path, landmarks: dict[str, list[float]]) -> Path:
    image = tmp_path / "sub-01_T1w.nii.gz"
    affine = np.array(
        [
            [-2.0, 0.0, 0.0, 10.0],
            [0.0, 3.0, 0.0, -20.0],
            [0.0, 0.0, -4.0, 30.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    nib.save(nib.Nifti1Image(np.zeros((2, 2, 2)), affine), image)
    image.with_name("sub-01_T1w.json").write_text(
        json.dumps({"AnatomicalLandmarkCoordinates": landmarks}),
        encoding="utf-8",
    )
    return image


def install_fake_3dcopy(monkeypatch: pytest.MonkeyPatch) -> list[list[str]]:
    calls: list[list[str]] = []

    def fake_run(args: list[str], *, check: bool) -> subprocess.CompletedProcess[str]:
        assert check is True
        calls.append(args)
        prefix = Path(args[2])
        Path(f"{prefix}+orig.HEAD").write_text(
            "type = string-attribute\n"
            "name = BYTEORDER_STRING\n"
            "count = 10\n"
            "'LSB_FIRST~\n",
            encoding="utf-8",
        )
        Path(f"{prefix}+orig.BRIK").write_bytes(b"brik")
        return subprocess.CompletedProcess(args, 0)

    monkeypatch.setattr(fiducials.shutil, "which", lambda command: "/afni/3dcopy")
    monkeypatch.setattr(fiducials.subprocess, "run", fake_run)
    return calls


def test_conversion_uses_the_input_affine_and_writes_afni_tags(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = write_bids_t1w(
        tmp_path,
        {"NAS": [1, 2, 3], "LPA": [4, 5, 6], "RPA": [7, 8, 9]},
    )
    calls = install_fake_3dcopy(monkeypatch)

    output_dir = tmp_path / "afni"
    brik, head = fiducials.convert_json_fids_to_head(image, output_dir)

    assert calls == [["/afni/3dcopy", str(image), str(output_dir / "sub-01_T1w")]]
    assert brik == output_dir / "sub-01_T1w+orig.BRIK"
    assert head == output_dir / "sub-01_T1w+orig.HEAD"
    header = head.read_text(encoding="utf-8")
    assert header.count("name = BYTEORDER_STRING") == 1
    assert "-8.000000\t14.000000\t18.000000\t0\t0" in header
    assert "-2.000000\t5.000000\t6.000000\t0\t0" in header
    assert "4.000000\t-4.000000\t-6.000000\t0\t0" in header
    assert "'Nasion~Left Ear~Right Ear~~~~~" in header


def test_existing_outputs_require_overwrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = write_bids_t1w(
        tmp_path,
        {"NAS": [0, 0, 0], "LPA": [0, 0, 0], "RPA": [0, 0, 0]},
    )
    calls = install_fake_3dcopy(monkeypatch)
    fiducials.convert_json_fids_to_head(image)

    with pytest.raises(
        fiducials.FiducialConversionError, match="output already exists"
    ):
        fiducials.convert_json_fids_to_head(image)
    assert len(calls) == 1

    fiducials.convert_json_fids_to_head(image, overwrite=True)
    assert len(calls) == 2


@pytest.mark.parametrize(
    "landmarks",
    [
        {"NAS": [1, 2, 3], "LPA": [4, 5, 6]},
        {"NAS": [1, 2], "LPA": [4, 5, 6], "RPA": [7, 8, 9]},
        {"NAS": [1, 2, 3], "LPA": [4, 5, 6], "RPA": [7, 8, float("nan")]},
    ],
)
def test_invalid_landmarks_are_rejected_before_running_afni(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    landmarks: dict[str, list[float]],
) -> None:
    image = write_bids_t1w(tmp_path, landmarks)
    monkeypatch.setattr(
        fiducials.subprocess,
        "run",
        lambda *args, **kwargs: pytest.fail("3dcopy should not run"),
    )

    with pytest.raises(fiducials.FiducialConversionError, match="NAS, LPA, and RPA"):
        fiducials.convert_json_fids_to_head(image)


def test_missing_afni_is_reported_at_conversion_time(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = write_bids_t1w(
        tmp_path,
        {"NAS": [0, 0, 0], "LPA": [0, 0, 0], "RPA": [0, 0, 0]},
    )
    monkeypatch.setattr(fiducials.shutil, "which", lambda command: None)

    with pytest.raises(fiducials.FiducialConversionError, match="3dcopy"):
        fiducials.convert_json_fids_to_head(image)


def test_input_filename_and_sidecar_are_validated(tmp_path: Path) -> None:
    with pytest.raises(fiducials.FiducialConversionError, match="end with .nii.gz"):
        fiducials.convert_json_fids_to_head(tmp_path / "sub-01_T1w.nii")

    image = tmp_path / "sub-01_T1w.nii.gz"
    nib.save(nib.Nifti1Image(np.zeros((2, 2, 2)), np.eye(4)), image)
    with pytest.raises(fiducials.FiducialConversionError, match="JSON sidecar"):
        fiducials.convert_json_fids_to_head(image)


def test_singular_affine_is_rejected_before_running_afni(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = write_bids_t1w(
        tmp_path,
        {"NAS": [0, 0, 0], "LPA": [0, 0, 0], "RPA": [0, 0, 0]},
    )
    monkeypatch.setattr(
        fiducials.nib,
        "load",
        lambda path: SimpleNamespace(affine=np.zeros((4, 4))),
    )
    monkeypatch.setattr(
        fiducials.subprocess,
        "run",
        lambda *args, **kwargs: pytest.fail("3dcopy should not run"),
    )

    with pytest.raises(fiducials.FiducialConversionError, match="nonsingular"):
        fiducials.convert_json_fids_to_head(image)


def test_help_does_not_require_afni(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(fiducials.shutil, "which", lambda command: None)
    monkeypatch.setattr("sys.argv", ["convert_json_fids_to_head", "--help"])

    with pytest.raises(SystemExit) as exit_info:
        fiducials.main()
    assert exit_info.value.code == 0
    assert "BIDS T1w NIfTI" in capsys.readouterr().out
