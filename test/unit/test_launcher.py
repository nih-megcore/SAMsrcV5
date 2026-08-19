from __future__ import annotations

import sys
from pathlib import Path

import pytest

from samsrcv5 import fiducials, launcher


pytestmark = pytest.mark.unit


def test_orthohull_argument_detection_preserves_optional_tag_file() -> None:
    image = Path("subject_T1w.nii.gz")
    argv = ["orthohull", "-i0", str(image), "fiducials.tag"]
    assert launcher._orthohull_nifti_argument(argv) == (2, image)


@pytest.mark.parametrize("entry_point", [launcher.orthohull, launcher.orthohull_python])
@pytest.mark.parametrize("suffix", [".nii", ".nii.gz"])
def test_orthohull_converts_nifti_and_removes_temporary_dataset(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    entry_point,
    suffix: str,
) -> None:
    image = tmp_path / f"subject_T1w{suffix}"
    image.touch()
    original_argv = ["orthohull", "-i", "0", "-x", "subject", str(image)]
    observed: dict[str, object] = {}

    def fake_convert(source: Path, output_dir: str):
        temporary = Path(output_dir)
        head = temporary / "subject_T1w+orig.HEAD"
        brik = temporary / "subject_T1w+orig.BRIK"
        head.write_text("tagged", encoding="utf-8")
        brik.write_bytes(b"brik")
        observed["source"] = source
        observed["temporary"] = temporary
        return brik, head

    def fake_legacy(name: str) -> None:
        observed["legacy"] = name
        observed["argv"] = sys.argv.copy()
        assert Path(sys.argv[-1]).read_text(encoding="utf-8") == "tagged"

    monkeypatch.setattr(fiducials, "convert_json_fids_to_head", fake_convert)
    monkeypatch.setattr(launcher, "_legacy", fake_legacy)
    monkeypatch.setattr(sys, "argv", original_argv)

    entry_point()

    assert observed["source"] == image
    assert observed["legacy"] == "orthohull.py"
    converted_argv = observed["argv"]
    assert converted_argv[:-1] == original_argv[:-1]
    assert converted_argv[-1].endswith("subject_T1w+orig.HEAD")
    assert sys.argv is original_argv
    assert not observed["temporary"].exists()


def test_orthohull_conversion_error_is_clear_and_skips_legacy_workflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = tmp_path / "subject_T1w.nii"
    image.touch()
    observed: dict[str, Path] = {}

    def fail_conversion(_source: Path, output_dir: str):
        observed["temporary"] = Path(output_dir)
        raise fiducials.FiducialConversionError("JSON sidecar does not exist")

    monkeypatch.setattr(fiducials, "convert_json_fids_to_head", fail_conversion)
    monkeypatch.setattr(
        launcher,
        "_legacy",
        lambda _name: pytest.fail("legacy workflow must not run"),
    )
    monkeypatch.setattr(sys, "argv", ["orthohull", str(image)])

    with pytest.raises(SystemExit, match="orthohull: error: JSON sidecar"):
        launcher.orthohull()

    assert not observed["temporary"].exists()


def test_orthohull_removes_temporary_dataset_when_legacy_workflow_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = tmp_path / "subject_T1w.nii.gz"
    image.touch()
    observed: dict[str, Path] = {}

    def fake_convert(_source: Path, output_dir: str):
        temporary = Path(output_dir)
        observed["temporary"] = temporary
        head = temporary / "subject_T1w+orig.HEAD"
        brik = temporary / "subject_T1w+orig.BRIK"
        head.touch()
        brik.touch()
        return brik, head

    monkeypatch.setattr(fiducials, "convert_json_fids_to_head", fake_convert)
    def fail_legacy(_name: str) -> None:
        raise RuntimeError("AFNI failed")

    monkeypatch.setattr(launcher, "_legacy", fail_legacy)
    monkeypatch.setattr(sys, "argv", ["orthohull", str(image)])

    with pytest.raises(RuntimeError, match="AFNI failed"):
        launcher.orthohull()

    assert not observed["temporary"].exists()


def test_orthohull_leaves_afni_input_unchanged(monkeypatch: pytest.MonkeyPatch) -> None:
    original_argv = ["orthohull", "-i0", "subject+orig.HEAD"]
    observed: dict[str, object] = {}

    monkeypatch.setattr(
        fiducials,
        "convert_json_fids_to_head",
        lambda *_args, **_kwargs: pytest.fail("AFNI input must not be converted"),
    )
    monkeypatch.setattr(
        launcher,
        "_legacy",
        lambda name: observed.update(name=name, argv=sys.argv.copy()),
    )
    monkeypatch.setattr(sys, "argv", original_argv)

    launcher.orthohull()

    assert observed == {"name": "orthohull.py", "argv": original_argv}
