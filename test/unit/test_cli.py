from __future__ import annotations

import math
import subprocess
from pathlib import Path

import pytest


pytestmark = pytest.mark.unit
ROOT = Path(__file__).resolve().parents[2]
BIN = ROOT / "bin"


def run(*args: str, input_text: str | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        input=input_text,
        text=True,
        capture_output=True,
        check=False,
    )


def test_expected_build_artifacts_exist() -> None:
    expected = {
        "1dstats",
        "OPMsim",
        "meshnorm",
        "orthohull",
        "sam_3d",
        "sam_4d",
        "sam_cov",
        "sam_ers",
        "sam_power",
        "sam_simulate",
        "sam_wts",
    }
    missing = sorted(name for name in expected if not (BIN / name).exists())
    assert not missing, f"missing installed programs: {missing}"


@pytest.mark.parametrize("program", ["sam_cov", "sam_wts", "sam_3d", "sam_4d", "sam_ers", "sam_power"])
def test_registered_program_help(program: str) -> None:
    result = run(str(BIN / program), "-h")
    assert result.returncode == 0
    output = result.stdout + result.stderr
    assert "Usage:" in output
    assert "Version 5.0" in output
    assert "-i_SAMdir SAMDIR" in output
    assert "-o_SAMdir SAMDIR" in output


def test_bad_numeric_argument_is_rejected() -> None:
    result = run(str(BIN / "sam_cov"), "--CovBand", "bad", "70")
    assert result.returncode != 0
    assert "badly formed number" in result.stderr


def test_sam_ers_uses_current_parameter_parser() -> None:
    result = run(str(BIN / "sam_ers"), "--TimeStep", "bad")
    assert result.returncode != 0
    assert "badly formed number" in result.stderr


@pytest.mark.parametrize("program", ["sam_3d", "sam_4d", "sam_ers", "sam_power"])
def test_image_metric_is_available_on_command_line(program: str) -> None:
    result = run(str(BIN / program), "--help")
    assert result.returncode == 0
    assert "--ImageMetric" in result.stdout + result.stderr


@pytest.mark.parametrize(
    ("program", "leading"),
    [
        ("sam_3d", ("--ImageMetric", "Power")),
        ("sam_4d", ("--Marker", "stim", "-0.1", "0.3", "TRUE")),
        ("sam_wts", ("--Model", "Nolte")),
        ("sam_wts", ("--ImageFormat", "ORIG")),
    ],
)
def test_variable_length_options_do_not_consume_following_options(
    program: str, leading: tuple[str, ...]
) -> None:
    result = run(str(BIN / program), *leading, "--CovBand", "bad", "70")
    assert result.returncode != 0
    assert "--CovBand: badly formed number 'bad'" in result.stderr


def test_complex_image_metric_does_not_consume_following_options() -> None:
    result = run(
        str(BIN / "sam_4d"),
        "--ImageMetric",
        "RankVectorEntropy",
        "0.01",
        "3",
        "--CovBand",
        "bad",
        "70",
    )
    assert result.returncode != 0
    assert "--CovBand: badly formed number 'bad'" in result.stderr


def test_marker_option_can_be_repeated() -> None:
    result = run(
        str(BIN / "sam_4d"),
        "--Marker",
        "stim",
        "-0.1",
        "0.3",
        "TRUE",
        "--Marker",
        "control",
        "-0.2",
        "0.4",
        "FALSE",
        "--TimeStep",
        "bad",
    )
    assert result.returncode != 0
    assert "--TimeStep: badly formed number 'bad'" in result.stderr


def test_sam_directory_option_requires_an_argument() -> None:
    result = run(str(BIN / "sam_cov"), "-i_SAMdir")
    assert result.returncode != 0
    assert "-i_SAMdir requires an argument" in result.stderr


def test_1dstats_known_sample() -> None:
    result = run(str(BIN / "1dstats"), "-q", input_text="1\n2\n3\n4\n")
    assert result.returncode == 0
    values = [float(value) for value in result.stdout.split()]
    assert values[0] == 4
    assert values[1] == pytest.approx(2.5)
    assert values[2] == pytest.approx(5.0 / 3.0)
    assert values[3:7] == pytest.approx([10.0, 1.0, 4.0, 1.0])
    assert values[7] == pytest.approx(math.sqrt(5.0 / 3.0))
    assert values[8] == pytest.approx(math.sqrt(5.0 / 3.0) / 2.0)
