from __future__ import annotations

import numpy as np
import pytest

mne = pytest.importorskip("mne")

from samsrcv5.mne import SAMNoise, estimate_sam_noise, make_sam_beamformer


pytestmark = pytest.mark.unit


def _info(channel_count: int):
    return mne.create_info(
        [f"MEG{index:03d}" for index in range(channel_count)],
        sfreq=100.0,
        ch_types="mag",
    )


def _covariance(info, data, *, projs=None, bads=None):
    return mne.Covariance(
        np.asarray(data, dtype=float),
        info["ch_names"],
        [] if bads is None else bads,
        [] if projs is None else projs,
        nfree=100,
    )


def _forward(info, gain, *, orientation_count: int, source_type: str = "vol"):
    gain = np.asarray(gain, dtype=float)
    source_count = gain.shape[1] // orientation_count
    if source_type == "surf":
        assert source_count == 2
        spaces = [
            {"type": "surf", "vertno": np.array([0]), "subject_his_id": "sample"},
            {"type": "surf", "vertno": np.array([1]), "subject_his_id": "sample"},
        ]
    else:
        spaces = [
            {
                "type": source_type,
                "vertno": np.arange(source_count),
                "subject_his_id": "sample",
            }
        ]
    source_spaces = mne.SourceSpaces(spaces)
    if orientation_count == 3:
        source_nn = np.tile(np.eye(3), (source_count, 1))
    else:
        source_nn = np.tile(np.array([[0.0, 0.0, 1.0]]), (source_count, 1))
    return mne.Forward(
        {
            "sol": {"data": gain, "row_names": list(info["ch_names"])},
            "nsource": source_count,
            "nchan": len(info["ch_names"]),
            "source_rr": np.column_stack(
                (np.arange(source_count), np.zeros((source_count, 2)))
            ),
            "source_nn": source_nn,
            "src": source_spaces,
            "info": info.copy(),
        }
    )


def test_fixed_orientation_matches_samsolve_and_applies_to_raw():
    info = _info(2)
    covariance = _covariance(info, np.diag([2.0, 4.0]))
    forward = _forward(info, [[1.0], [1.0]], orientation_count=1)

    filters = make_sam_beamformer(info, forward, covariance)

    assert isinstance(filters, mne.beamformer.Beamformer)
    assert filters["kind"] == "SAM"
    np.testing.assert_allclose(filters["weights"], [[2.0 / 3.0, 1.0 / 3.0]])
    np.testing.assert_allclose(filters["sam_condition_numbers"], [0.0])
    assert filters["src_type"] == "volume"

    data = np.arange(20.0).reshape(2, 10)
    raw = mne.io.RawArray(data, info, verbose=False)
    source = mne.beamformer.apply_lcmv_raw(raw, filters, verbose=False)
    assert isinstance(source, mne.VolSourceEstimate)
    np.testing.assert_allclose(source.data, filters["weights"] @ data)


def test_free_orientation_uses_separate_orientation_covariance():
    info = _info(3)
    data_cov = _covariance(info, np.diag([2.0, 3.0, 4.0]))
    orient_cov = _covariance(info, np.diag([1.0, 2.0, 4.0]))
    forward = _forward(info, np.eye(3), orientation_count=3)

    filters = make_sam_beamformer(
        info,
        forward,
        data_cov,
        orient_cov=orient_cov,
    )

    assert filters["pick_ori"] == "max-power"
    np.testing.assert_allclose(filters["max_power_ori"], [[0.0, 0.0, 1.0]])
    np.testing.assert_allclose(filters["sam_orientations"], [[0.0, 0.0, 1.0]])
    np.testing.assert_allclose(filters["sam_condition_numbers"], [4.0])
    np.testing.assert_allclose(filters["weights"], [[0.0, 0.0, 1.0]])


def test_free_orientation_reduces_silent_orientation_rank():
    info = _info(3)
    covariance = _covariance(info, np.diag([1.0, 2.0, 3.0]))
    forward = _forward(
        info,
        np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]]),
        orientation_count=3,
    )

    filters = make_sam_beamformer(info, forward, covariance)

    np.testing.assert_allclose(filters["max_power_ori"], [[0.0, 1.0, 0.0]])
    np.testing.assert_allclose(filters["sam_condition_numbers"], [2.0])
    np.testing.assert_allclose(filters["weights"], [[0.0, 1.0, 0.0]])


def test_noise_estimation_matches_sam_spectrum_search():
    info = _info(6)
    covariance = _covariance(info, np.diag([10.0, 8.0, 6.0, 4.0, 2.0, 1.0]))

    noise = estimate_sam_noise(covariance, 20.0)

    assert noise == SAMNoise(variance=4.0, bandwidth_hz=20.0)


def test_proportional_regularization_and_noise_normalization():
    info = _info(2)
    covariance = _covariance(info, np.diag([2.0, 4.0]))
    forward = _forward(info, [[1.0], [1.0]], orientation_count=1)

    filters = make_sam_beamformer(
        info,
        forward,
        covariance,
        data_noise=SAMNoise(1.0),
        mu=1.0,
        mu_kind="proportional",
        normalize=True,
    )

    expected = np.array([5.0 / 8.0, 3.0 / 8.0]) / np.sqrt(2.0)
    np.testing.assert_allclose(filters["weights"], [expected])
    assert filters["sam_effective_noise"] == pytest.approx(2.0)
    assert filters["sam_data_loading"] == pytest.approx(1.0)


def test_projection_is_preserved_for_mne_application():
    info = _info(3)
    projection = mne.Projection(
        data={
            "col_names": list(info["ch_names"]),
            "row_names": None,
            "data": np.array([[1.0, 0.0, 0.0]]),
            "nrow": 1,
            "ncol": 3,
        },
        desc="remove first channel",
        active=False,
    )
    template = mne.io.RawArray(np.zeros((3, 2)), info, verbose=False)
    template.add_proj([projection], verbose=False)
    info = template.info
    covariance = _covariance(info, np.eye(3), projs=info["projs"])
    forward = _forward(info, [[0.0], [1.0], [1.0]], orientation_count=1)

    filters = make_sam_beamformer(info, forward, covariance)
    data = np.arange(30.0).reshape(3, 10)
    raw = mne.io.RawArray(data, info, verbose=False)
    source = mne.beamformer.apply_lcmv_raw(raw, filters, verbose=False)

    assert filters["is_ssp"] is True
    np.testing.assert_allclose(source.data, filters["weights"] @ filters["proj"] @ data)


def test_surface_source_metadata_is_retained():
    info = _info(2)
    covariance = _covariance(info, np.eye(2))
    forward = _forward(info, np.eye(2), orientation_count=1, source_type="surf")

    filters = make_sam_beamformer(info, forward, covariance)

    assert filters["src_type"] == "surface"
    assert filters["subject"] == "sample"
    assert [vertex.tolist() for vertex in filters["vertices"]] == [[0], [1]]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"normalize": True}, "positive data_noise"),
        ({"mu": 1.0}, "data_noise is required"),
        (
            {"data_noise": SAMNoise(1.0), "mu": 1.0},
            "bandwidth_hz is required",
        ),
        ({"n_nulls": 2}, "smaller than the channel count"),
    ],
)
def test_invalid_noise_and_rank_options_raise(kwargs, message):
    info = _info(2)
    covariance = _covariance(info, np.eye(2))
    forward = _forward(info, [[1.0], [1.0]], orientation_count=1)

    with pytest.raises(ValueError, match=message):
        make_sam_beamformer(info, forward, covariance, **kwargs)


def test_noise_metadata_requires_sam_noise():
    info = _info(2)
    covariance = _covariance(info, np.eye(2))
    forward = _forward(info, [[1.0], [1.0]], orientation_count=1)

    with pytest.raises(TypeError, match="data_noise must be a SAMNoise"):
        make_sam_beamformer(info, forward, covariance, data_noise=1.0)


def test_mismatched_channels_and_projections_raise():
    info = _info(2)
    incomplete_covariance = mne.Covariance(
        np.array([[1.0]]), [info["ch_names"][0]], [], [], nfree=10
    )
    forward = _forward(info, [[1.0], [1.0]], orientation_count=1)
    with pytest.raises(ValueError, match="missing channels"):
        make_sam_beamformer(info, forward, incomplete_covariance)

    projection = mne.Projection(
        data={
            "col_names": list(info["ch_names"]),
            "row_names": None,
            "data": np.array([[1.0, 0.0]]),
            "nrow": 1,
            "ncol": 2,
        },
        desc="mismatch",
        active=False,
    )
    projected_covariance = _covariance(info, np.eye(2), projs=[projection])
    with pytest.raises(ValueError, match="projections do not match"):
        make_sam_beamformer(info, forward, projected_covariance)


def test_construction_performs_no_file_io(tmp_path, monkeypatch):
    info = _info(2)
    covariance = _covariance(info, np.eye(2))
    forward = _forward(info, [[1.0], [1.0]], orientation_count=1)
    monkeypatch.chdir(tmp_path)

    make_sam_beamformer(info, forward, covariance)

    assert list(tmp_path.iterdir()) == []
