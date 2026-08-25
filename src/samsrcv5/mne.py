"""In-memory SAM beamformer construction for MNE-Python objects.

This module intentionally does not call the native SAM executables or any of
the SAM image writers.  It implements the small linear-algebra boundary needed
to turn an MNE forward solution and covariance matrices into an MNE
``Beamformer`` object.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import linalg

try:
    import mne
    from mne.beamformer import Beamformer
    from mne.io.constants import FIFF
    from mne.proj import make_projector
except ModuleNotFoundError as error:  # pragma: no cover - exercised without extra
    raise ModuleNotFoundError(
        "The samsrcv5 MNE interface requires the optional dependency; "
        "install it with 'samsrcv5[mne]'."
    ) from error


MuKind = Literal["additive", "proportional"]


@dataclass(frozen=True)
class SAMNoise:
    """Noise metadata used by legacy SAM regularization.

    Parameters
    ----------
    variance
        Sensor noise variance in the squared units of the covariance matrix.
    bandwidth_hz
        Effective filter bandwidth.  This is required for additive ``mu``,
        whose public value is expressed in fT/sqrt(Hz), as in ``sam_wts``.
    """

    variance: float
    bandwidth_hz: float | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.variance) or self.variance < 0.0:
            raise ValueError("noise variance must be finite and non-negative")
        if self.bandwidth_hz is not None and (
            not np.isfinite(self.bandwidth_hz) or self.bandwidth_hz <= 0.0
        ):
            raise ValueError("noise bandwidth_hz must be finite and positive")


def _covariance_array(covariance: mne.Covariance, name: str) -> np.ndarray:
    names = list(covariance.ch_names)
    data = np.asarray(covariance["data"], dtype=float)
    if bool(covariance.get("diag", False)) or data.ndim == 1:
        if data.shape != (len(names),):
            raise ValueError(f"{name} diagonal does not match its channel names")
        data = np.diag(data)
    if data.shape != (len(names), len(names)):
        raise ValueError(f"{name} must be a square channel-by-channel matrix")
    if not np.all(np.isfinite(data)):
        raise ValueError(f"{name} contains non-finite values")
    scale = float(np.max(np.abs(data))) if data.size else 0.0
    tolerance = 100.0 * np.finfo(float).eps * scale
    if float(np.max(np.abs(data - data.T))) > tolerance:
        raise ValueError(f"{name} must be symmetric")
    data = 0.5 * (data + data.T)
    eigenvalues = np.linalg.eigvalsh(data)
    if eigenvalues.size and eigenvalues[0] < -max(tolerance, scale * 1.0e-12):
        raise ValueError(f"{name} must be positive semidefinite")
    return data


def estimate_sam_noise(
    covariance: mne.Covariance,
    bandwidth_hz: float,
    *,
    edge_exclusion: int = 3,
) -> SAMNoise:
    """Estimate the SAM sensor-noise floor from a covariance spectrum.

    This follows the singular-value derivative search used by ``sam_cov``.
    It does not modify the covariance object.
    """

    if not isinstance(covariance, mne.Covariance):
        raise TypeError("covariance must be an mne.Covariance")
    if not np.isfinite(bandwidth_hz) or bandwidth_hz <= 0.0:
        raise ValueError("bandwidth_hz must be finite and positive")
    if not isinstance(edge_exclusion, int) or edge_exclusion < 3:
        raise ValueError("edge_exclusion must be an integer of at least 3")

    values = np.linalg.svd(
        _covariance_array(covariance, "covariance"), compute_uv=False
    )
    count = len(values)
    if count < 5 or edge_exclusion > count - 2:
        raise ValueError(
            "covariance has too few channels for the requested edge exclusion"
        )

    best_index: int | None = None
    best_delta = np.inf
    for index in range(count - edge_exclusion, 1, -1):
        delta = values[index - 2] - values[index + 2]
        if delta < best_delta:
            best_delta = float(delta)
            best_index = index
    if best_index is None:  # Defensive; validated dimensions make this unreachable.
        raise ValueError("could not estimate a noise floor from the covariance")
    return SAMNoise(float(values[best_index]), float(bandwidth_hz))


def _sam_pinv(
    matrix: np.ndarray,
    *,
    rcond: float,
    n_nulls: int,
) -> tuple[np.ndarray, int]:
    u, singular_values, vh = np.linalg.svd(matrix, full_matrices=False)
    keep = singular_values > singular_values[0] * rcond
    if n_nulls:
        keep[-n_nulls:] = False
    inverse_values = np.zeros_like(singular_values)
    inverse_values[keep] = 1.0 / singular_values[keep]
    inverse = (vh.T * inverse_values) @ u.T
    return 0.5 * (inverse + inverse.T), int(np.count_nonzero(keep))


def _regularize(
    covariance: np.ndarray,
    noise: SAMNoise | None,
    *,
    mu: float | None,
    mu_kind: MuKind,
    tesla_units: bool,
    name: str,
) -> tuple[np.ndarray, float | None, float]:
    effective_noise = None if noise is None else float(noise.variance)
    loading = 0.0
    if mu is not None:
        if noise is None:
            raise ValueError(f"{name}_noise is required when mu is specified")
        if mu_kind == "additive":
            if not tesla_units:
                raise ValueError(
                    "additive mu in fT/sqrt(Hz) requires MEG channels in tesla"
                )
            if noise.bandwidth_hz is None:
                raise ValueError(
                    f"{name}_noise.bandwidth_hz is required for additive mu"
                )
            loading = noise.bandwidth_hz * (mu * 1.0e-15) ** 2
        else:
            loading = noise.variance * mu
        effective_noise += loading

    result = covariance.copy()
    if loading:
        result.flat[:: result.shape[0] + 1] += loading
    return result, effective_noise, float(loading)


def _projection_matrix(projections: list, channel_names: list[str]) -> np.ndarray:
    projection, _, _ = make_projector(projections, channel_names)
    return np.asarray(projection, dtype=float)


def _require_matching_projection(
    expected: np.ndarray,
    projections: list,
    channel_names: list[str],
    name: str,
) -> None:
    actual = _projection_matrix(projections, channel_names)
    if not np.allclose(
        actual,
        expected,
        rtol=1.0e-12,
        atol=np.finfo(float).eps,
    ):
        raise ValueError(f"{name} projections do not match info projections")


def _ordered_covariance(
    covariance: mne.Covariance,
    channel_names: list[str],
    name: str,
) -> np.ndarray:
    all_names = list(covariance.ch_names)
    missing = [channel for channel in channel_names if channel not in all_names]
    if missing:
        raise ValueError(f"{name} is missing channels: {', '.join(missing)}")
    indices = [all_names.index(channel) for channel in channel_names]
    full = _covariance_array(covariance, name)
    return full[np.ix_(indices, indices)]


def _source_metadata(
    forward: mne.Forward,
    source_count: int,
) -> tuple[list[np.ndarray], str, str | None]:
    source_spaces = forward["src"]
    vertices = [
        np.asarray(space["vertno"], dtype=int).copy() for space in source_spaces
    ]
    if sum(len(vertex) for vertex in vertices) != source_count:
        raise ValueError("forward source-space vertices do not match nsource")
    source_type = source_spaces.kind
    subject = getattr(source_spaces, "_subject", None)
    return vertices, source_type, subject


def _optimal_orientation(
    leadfield: np.ndarray,
    covariance_inverse: np.ndarray,
    basis: np.ndarray,
    *,
    rcond: float,
    source_index: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    projected = covariance_inverse @ leadfield
    source_power = leadfield.T @ projected
    noise_power = projected.T @ projected
    source_power = 0.5 * (source_power + source_power.T)
    noise_power = 0.5 * (noise_power + noise_power.T)

    noise_eigenvalues, noise_eigenvectors = np.linalg.eigh(noise_power)
    maximum = float(np.max(noise_eigenvalues))
    if maximum <= 0.0:
        raise np.linalg.LinAlgError(
            f"source {source_index} has a degenerate orientation leadfield"
        )
    keep = noise_eigenvalues > maximum * rcond
    if not np.any(keep):
        raise np.linalg.LinAlgError(
            f"source {source_index} has no usable orientation dimensions"
        )

    whitening = noise_eigenvectors[:, keep] / np.sqrt(noise_eigenvalues[keep])
    reduced = whitening.T @ source_power @ whitening
    reduced = 0.5 * (reduced + reduced.T)
    eigenvalues, eigenvectors = linalg.eigh(reduced, check_finite=False)
    local_orientation = whitening @ eigenvectors[:, -1]
    norm = float(np.linalg.norm(local_orientation))
    if not np.isfinite(norm) or norm == 0.0:
        raise np.linalg.LinAlgError(
            f"source {source_index} produced an invalid orientation"
        )
    local_orientation /= norm
    physical_orientation = local_orientation @ basis

    reference_normal = basis[2]
    alignment = float(np.dot(physical_orientation, reference_normal))
    if abs(alignment) <= np.finfo(float).eps:
        sign_index = int(np.argmax(np.abs(local_orientation)))
        alignment = float(local_orientation[sign_index])
    if alignment < 0.0:
        local_orientation *= -1.0
        physical_orientation *= -1.0

    positive = eigenvalues[eigenvalues > abs(eigenvalues[-1]) * rcond]
    condition = np.inf if len(positive) < 2 else float(positive[-1] / positive[0])
    return local_orientation, physical_orientation, condition


def make_sam_beamformer(
    info: mne.Info,
    forward: mne.Forward,
    data_cov: mne.Covariance,
    *,
    orient_cov: mne.Covariance | None = None,
    data_noise: SAMNoise | None = None,
    orient_noise: SAMNoise | None = None,
    mu: float | None = None,
    mu_kind: MuKind = "additive",
    normalize: bool = False,
    n_nulls: int = 0,
    rcond: float = 1.0e-15,
) -> Beamformer:
    """Construct native-style SAM weights entirely in memory.

    The input forward solution supplies the source grid and leadfield.  A
    separate orientation covariance can be supplied to reproduce SAM's
    ``Orient.cov`` behavior.  The returned scalar beamformer can be passed to
    MNE's ``apply_lcmv_*`` functions.
    """

    if not isinstance(info, mne.Info):
        raise TypeError("info must be an mne.Info")
    if not isinstance(forward, mne.Forward):
        raise TypeError("forward must be an mne.Forward")
    if not isinstance(data_cov, mne.Covariance):
        raise TypeError("data_cov must be an mne.Covariance")
    if orient_cov is not None and not isinstance(orient_cov, mne.Covariance):
        raise TypeError("orient_cov must be an mne.Covariance or None")
    if data_noise is not None and not isinstance(data_noise, SAMNoise):
        raise TypeError("data_noise must be a SAMNoise or None")
    if orient_noise is not None and not isinstance(orient_noise, SAMNoise):
        raise TypeError("orient_noise must be a SAMNoise or None")
    if mu_kind not in ("additive", "proportional"):
        raise ValueError("mu_kind must be 'additive' or 'proportional'")
    if mu is not None and (not np.isfinite(mu) or mu < 0.0):
        raise ValueError("mu must be finite and non-negative")
    if not isinstance(normalize, bool):
        raise TypeError("normalize must be a bool")
    if not isinstance(n_nulls, int) or n_nulls < 0:
        raise ValueError("n_nulls must be a non-negative integer")
    if not np.isfinite(rcond) or rcond <= 0.0 or rcond >= 1.0:
        raise ValueError("rcond must be finite and between zero and one")

    if orient_cov is None:
        orient_cov = data_cov
        if orient_noise is None:
            orient_noise = data_noise

    bads = set(info["bads"])
    bads.update(data_cov.get("bads", []))
    bads.update(orient_cov.get("bads", []))
    forward_info = forward.get("info")
    if forward_info is not None:
        bads.update(forward_info.get("bads", []))
    picks = mne.pick_types(
        info,
        meg=True,
        ref_meg=False,
        eeg=False,
        exclude=sorted(bads),
    )
    if len(picks) == 0:
        raise ValueError("info contains no usable non-reference MEG channels")
    channel_names = [info["ch_names"][pick] for pick in picks]
    channel_units = {info["chs"][pick]["unit"] for pick in picks}
    if len(channel_units) != 1:
        raise ValueError("SAM requires homogeneous MEG channel units")
    tesla_units = channel_units == {FIFF.FIFF_UNIT_T}

    forward_rows = list(forward["sol"].get("row_names") or [])
    if not forward_rows and forward_info is not None:
        forward_rows = list(forward_info["ch_names"])
    missing = [channel for channel in channel_names if channel not in forward_rows]
    if missing:
        raise ValueError(f"forward is missing channels: {', '.join(missing)}")
    if len(set(forward_rows)) != len(forward_rows):
        raise ValueError("forward channel names must be unique")
    forward_indices = [forward_rows.index(channel) for channel in channel_names]
    gain = np.asarray(forward["sol"]["data"], dtype=float)[forward_indices]
    if not np.all(np.isfinite(gain)):
        raise ValueError("forward leadfield contains non-finite values")

    data_matrix = _ordered_covariance(data_cov, channel_names, "data_cov")
    orient_matrix = _ordered_covariance(orient_cov, channel_names, "orient_cov")
    projector = _projection_matrix(info["projs"], channel_names)
    _require_matching_projection(
        projector, data_cov.get("projs", []), channel_names, "data_cov"
    )
    _require_matching_projection(
        projector, orient_cov.get("projs", []), channel_names, "orient_cov"
    )
    if forward_info is not None:
        _require_matching_projection(
            projector,
            forward_info.get("projs", []),
            channel_names,
            "forward",
        )
    gain = projector @ gain
    data_matrix = projector @ data_matrix @ projector.T
    orient_matrix = projector @ orient_matrix @ projector.T

    data_matrix, effective_noise, data_loading = _regularize(
        data_matrix,
        data_noise,
        mu=mu,
        mu_kind=mu_kind,
        tesla_units=tesla_units,
        name="data",
    )
    orient_matrix, effective_orient_noise, orient_loading = _regularize(
        orient_matrix,
        orient_noise,
        mu=mu,
        mu_kind=mu_kind,
        tesla_units=tesla_units,
        name="orient",
    )
    if normalize and (effective_noise is None or effective_noise <= 0.0):
        raise ValueError("normalize=True requires positive data_noise")

    channel_count = len(channel_names)
    if n_nulls >= channel_count:
        raise ValueError("n_nulls must be smaller than the channel count")
    data_inverse, data_rank = _sam_pinv(
        data_matrix, rcond=rcond, n_nulls=n_nulls
    )
    orient_inverse, orient_rank = _sam_pinv(
        orient_matrix, rcond=rcond, n_nulls=n_nulls
    )

    source_count = int(forward["nsource"])
    if source_count <= 0 or gain.shape[1] % source_count:
        raise ValueError("forward gain dimensions do not match nsource")
    orientation_count = gain.shape[1] // source_count
    if orientation_count not in (1, 3):
        raise ValueError("forward must have one or three orientations per source")
    source_rr = np.asarray(forward["source_rr"], dtype=float)
    if source_rr.shape != (source_count, 3) or not np.all(np.isfinite(source_rr)):
        raise ValueError("forward source_rr must have shape (nsource, 3)")
    source_nn = np.asarray(forward["source_nn"], dtype=float)
    if not np.all(np.isfinite(source_nn)):
        raise ValueError("forward source_nn must contain only finite values")

    weights = np.empty((source_count, channel_count), dtype=float)
    local_orientations: np.ndarray | None
    physical_orientations = np.empty((source_count, 3), dtype=float)
    condition_numbers = np.zeros(source_count, dtype=float)
    if orientation_count == 3:
        if source_nn.shape != (3 * source_count, 3):
            raise ValueError(
                "free-orientation forward source_nn must have shape (3*nsource, 3)"
            )
        bases = source_nn.reshape(source_count, 3, 3)
        local_orientations = np.empty((source_count, 3), dtype=float)
    else:
        if source_nn.shape == (3 * source_count, 3):
            fixed_normals = source_nn.reshape(source_count, 3, 3)[:, 2]
        elif source_nn.shape == (source_count, 3):
            fixed_normals = source_nn
        else:
            raise ValueError(
                "fixed-orientation forward source_nn must have shape (nsource, 3)"
            )
        norms = np.linalg.norm(fixed_normals, axis=1)
        if np.any(norms == 0.0) or not np.all(np.isfinite(norms)):
            raise ValueError(
                "fixed-orientation source normals must be finite and nonzero"
            )
        physical_orientations[:] = fixed_normals / norms[:, np.newaxis]
        local_orientations = None

    normalization = 1.0 if not normalize else 1.0 / np.sqrt(effective_noise)
    for source_index in range(source_count):
        start = source_index * orientation_count
        leadfield = gain[:, start : start + orientation_count]
        if orientation_count == 3:
            local, physical, condition = _optimal_orientation(
                leadfield,
                orient_inverse,
                bases[source_index],
                rcond=rcond,
                source_index=source_index,
            )
            local_orientations[source_index] = local
            physical_orientations[source_index] = physical
            condition_numbers[source_index] = condition
            forward_vector = leadfield @ local
        else:
            forward_vector = leadfield[:, 0]

        unnormalized = data_inverse @ forward_vector
        denominator = float(forward_vector @ unnormalized)
        denominator_scale = np.linalg.norm(forward_vector) * np.linalg.norm(
            unnormalized
        )
        if not np.isfinite(denominator) or abs(denominator) <= (
            np.finfo(float).eps * denominator_scale
        ):
            raise np.linalg.LinAlgError(
                f"source {source_index} has zero unit-gain denominator"
            )
        weights[source_index] = normalization * unnormalized / denominator

    vertices, source_type, subject = _source_metadata(forward, source_count)
    picked_data_cov = mne.pick_channels_cov(
        data_cov,
        include=channel_names,
        exclude=[],
        ordered=True,
        copy=True,
    )
    return Beamformer(
        kind="SAM",
        weights=weights,
        data_cov=picked_data_cov,
        noise_cov=None,
        whitener=None,
        weight_norm=None,
        pick_ori="max-power" if orientation_count == 3 else None,
        ch_names=channel_names,
        proj=projector,
        is_ssp=bool(info["projs"]),
        vertices=vertices,
        is_free_ori=False,
        n_sources=source_count,
        src_type=source_type,
        source_nn=source_nn.copy(),
        subject=subject,
        rank=data_rank,
        max_power_ori=None if local_orientations is None else local_orientations,
        inversion="single",
        sam_source_rr=source_rr.copy(),
        sam_orientations=physical_orientations,
        sam_condition_numbers=condition_numbers,
        sam_effective_noise=effective_noise,
        sam_effective_orient_noise=effective_orient_noise,
        sam_data_loading=data_loading,
        sam_orient_loading=orient_loading,
        sam_mu=mu,
        sam_mu_kind=mu_kind if mu is not None else None,
        sam_n_nulls=n_nulls,
        sam_rcond=float(rcond),
        sam_orient_rank=orient_rank,
    )


__all__ = ["SAMNoise", "estimate_sam_noise", "make_sam_beamformer"]
