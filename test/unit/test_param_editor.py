from __future__ import annotations

from pathlib import Path

import pytest

from samsrcv5.param_editor import (
    SPEC_BY_KEY,
    ParameterDocument,
    canonical_key,
    changed_updates,
    ensure_save_root,
    managed_parameter_path,
    serialize_new,
    validate_values,
)


pytestmark = pytest.mark.unit


def test_catalog_uses_current_names_and_program_profiles() -> None:
    assert canonical_key("marker12") == "Marker"
    assert canonical_key("inputsamd") == "InputSAMDirectory"
    assert canonical_key("Image") is None
    assert "sam_ers" in SPEC_BY_KEY["SmoothBand"].programs
    assert "sam_wts" in SPEC_BY_KEY["Model"].programs
    assert "PropMu" not in SPEC_BY_KEY
    assert "Surface" not in SPEC_BY_KEY


def test_document_round_trip_preserves_raw_content() -> None:
    source = (
        "# investigator comment\n"
        "CovBand 5 70  # acquisition-safe\n"
        "SiteSpecific keep this exactly\n"
        "\n"
        "Marker1 stim -0.1 0.3 TRUE\n\n"
    )
    document = ParameterDocument.parse(source)
    assert document.values() == {
        "CovBand": ["5 70"],
        "Marker": ["stim -0.1 0.3 TRUE"],
    }
    assert document.render({}) == source

    current = {"CovBand": ["1 40"], "Marker": ["stim -0.1 0.3 TRUE"]}
    rendered = document.render(changed_updates(document.values(), current))
    assert "CovBand 1 40  # acquisition-safe" in rendered
    assert "SiteSpecific keep this exactly" in rendered
    assert "Marker1 stim -0.1 0.3 TRUE" in rendered


def test_repeat_parameters_and_flags_serialize_canonically() -> None:
    text = serialize_new(
        {
            "Verbose": [""],
            "Marker": ["stim -0.1 0 TRUE", "stim 0 0.3 FALSE post"],
            "CovBand": ["5 70"],
        }
    )
    assert text == (
        "Verbose\n"
        "Marker stim -0.1 0 TRUE\n"
        "Marker stim 0 0.3 FALSE post\n"
        "CovBand 5 70\n"
    )


def test_validation_reports_syntax_errors_and_profile_warnings() -> None:
    values = {
        "Marker": ["stim -0.1 0.3 TRUE"],
        "CovType": ["ALL"],
        "CovBand": ["5 40"],
        "ImageBand": ["1 70"],
        "SmoothBand": ["0 20"],
        "TimeStep": ["0.01"],
        "ImageMetric": ["Signal"],
    }
    errors, warnings = validate_values(values, "sam_ers")
    assert "CovType: sam_ers supports only GLOBAL or SUM" in errors
    assert "ImageBand must lie within CovBand" in errors
    assert not warnings

    errors, warnings = validate_values({"CovBand": ["5 70"]}, "sam_wts")
    assert not errors
    assert any("Model" in warning for warning in warnings)
    assert any("XYZ grid" in warning for warning in warnings)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("Model", "SingleSphere 0 0 4"),
        ("Model", "Nolte 16"),
        ("Mu", "*1.5"),
        ("ImageFormat", "TLRC 2"),
        ("PrefixLength", "_"),
    ],
)
def test_complex_values_validate(key: str, value: str) -> None:
    errors, _warnings = validate_values({key: [value]}, "sam_cov")
    assert not errors


def test_sam_3d_accepts_case_insensitive_power_metric() -> None:
    errors, _warnings = validate_values({"ImageMetric": ["power"]}, "sam_3d")
    assert not errors


def test_managed_parameter_path_is_flat_and_adds_extension(tmp_path: Path) -> None:
    assert managed_parameter_path(tmp_path, "analysis") == tmp_path / "analysis.param"
    assert managed_parameter_path(tmp_path, "../elsewhere.param") == tmp_path / "elsewhere.param"
    with pytest.raises(ValueError):
        managed_parameter_path(tmp_path, "")


def test_save_root_is_created_under_home(tmp_path: Path) -> None:
    assert ensure_save_root(tmp_path) == tmp_path / "samparams"
    assert (tmp_path / "samparams").is_dir()


def test_gui_module_import_does_not_require_a_display() -> None:
    from samsrcv5 import param_gui

    assert callable(param_gui.main)
