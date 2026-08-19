"""Convert BIDS T1w fiducials to an AFNI BRIK/HEAD dataset."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np


LANDMARK_NAMES = ("NAS", "LPA", "RPA")
AFNI_TAG_NAMES = ("TAGSET_NUM", "TAGSET_FLOATS", "TAGSET_LABELS")

__all__ = ["FiducialConversionError", "convert_json_fids_to_head"]


class FiducialConversionError(RuntimeError):
    """Raised when a BIDS T1w image cannot be converted safely."""


def _input_paths(t1w_nii: str | Path) -> tuple[Path, Path, str]:
    image = Path(t1w_nii)
    if not image.name.endswith(".nii.gz"):
        raise FiducialConversionError("input filename must end with .nii.gz")
    stem = image.name.removesuffix(".nii.gz")
    sidecar = image.with_name(f"{stem}.json")
    if not image.is_file():
        raise FiducialConversionError(f"NIfTI file does not exist: {image}")
    if not sidecar.is_file():
        raise FiducialConversionError(f"JSON sidecar does not exist: {sidecar}")
    return image, sidecar, stem


def _load_landmarks(sidecar: Path) -> np.ndarray:
    try:
        metadata = json.loads(sidecar.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise FiducialConversionError(
            f"could not read JSON sidecar {sidecar}: {error}"
        ) from error

    try:
        landmarks = metadata["AnatomicalLandmarkCoordinates"]
        coordinates = np.asarray(
            [landmarks[name] for name in LANDMARK_NAMES], dtype=float
        )
    except (KeyError, TypeError, ValueError) as error:
        raise FiducialConversionError(
            "AnatomicalLandmarkCoordinates must contain numeric "
            "NAS, LPA, and RPA values"
        ) from error

    if coordinates.shape != (3, 3) or not np.isfinite(coordinates).all():
        raise FiducialConversionError(
            "NAS, LPA, and RPA must each contain three finite voxel coordinates"
        )
    return coordinates


def _load_affine(image: Path) -> np.ndarray:
    try:
        affine = np.asarray(nib.load(image).affine, dtype=float)
    except Exception as error:
        raise FiducialConversionError(
            f"could not read NIfTI image {image}: {error}"
        ) from error

    if affine.shape != (4, 4) or not np.isfinite(affine).all():
        raise FiducialConversionError("NIfTI affine must be a finite 4-by-4 matrix")
    if np.isclose(np.linalg.det(affine[:3, :3]), 0.0):
        raise FiducialConversionError(
            "NIfTI affine must have a nonsingular spatial transform"
        )
    return affine


def _landmarks_in_afni_lps(
    affine: np.ndarray, voxel_coordinates: np.ndarray
) -> np.ndarray:
    ras_coordinates = nib.affines.apply_affine(affine, voxel_coordinates)
    return np.round(ras_coordinates * np.array([-1.0, -1.0, 1.0]), decimals=6)


def _format_landmark(coordinate: np.ndarray) -> str:
    values = [*(f"{value:.6f}" for value in coordinate), "0", "0"]
    return "\t".join(values)


def _tagset_text(coordinates: np.ndarray) -> str:
    nas, lpa, rpa = (_format_landmark(coordinate) for coordinate in coordinates)
    return (
        "\n"
        "type = integer-attribute\n"
        "name = TAGSET_NUM\n"
        "count = 2\n"
        " 3 5\n"
        "\n"
        "type = float-attribute\n"
        "name = TAGSET_FLOATS\n"
        "count = 15\n"
        f"{nas}\n"
        f"{lpa}\n"
        f"{rpa}\n"
        "\n"
        "type = string-attribute\n"
        "name = TAGSET_LABELS\n"
        "count = 30\n"
        "'Nasion~Left Ear~Right Ear~~~~~\n"
    )


def _append_tagset(head: Path, coordinates: np.ndarray) -> None:
    try:
        header = head.read_text(encoding="utf-8")
    except OSError as error:
        raise FiducialConversionError(
            f"could not read AFNI header {head}: {error}"
        ) from error

    for name in AFNI_TAG_NAMES:
        if re.search(rf"^name\s*=\s*{re.escape(name)}\s*$", header, re.MULTILINE):
            raise FiducialConversionError(
                f"AFNI header already contains {name}: {head}"
            )

    try:
        with head.open("a", encoding="utf-8") as stream:
            if header and not header.endswith("\n"):
                stream.write("\n")
            stream.write(_tagset_text(coordinates))
    except OSError as error:
        raise FiducialConversionError(
            f"could not update AFNI header {head}: {error}"
        ) from error


def convert_json_fids_to_head(
    t1w_nii: str | Path,
    output_dir: str | Path | None = None,
    *,
    overwrite: bool = False,
) -> tuple[Path, Path]:
    """Create an AFNI dataset carrying fiducials from a BIDS JSON sidecar.

    Landmark coordinates are interpreted as voxel coordinates for ``t1w_nii``.
    They are transformed by that image's affine and converted from RAS world
    coordinates to the LPS coordinates stored by AFNI.

    Returns the generated ``(BRIK, HEAD)`` paths.
    """

    image, sidecar, stem = _input_paths(t1w_nii)
    voxel_coordinates = _load_landmarks(sidecar)
    affine = _load_affine(image)
    lps_coordinates = _landmarks_in_afni_lps(affine, voxel_coordinates)

    destination = Path(output_dir) if output_dir is not None else image.parent
    prefix = destination / stem
    head = Path(f"{prefix}+orig.HEAD")
    brik = Path(f"{prefix}+orig.BRIK")
    compressed_brik = Path(f"{brik}.gz")
    outputs = (head, brik, compressed_brik)

    existing = [path for path in outputs if path.exists()]
    if existing and not overwrite:
        names = ", ".join(str(path) for path in existing)
        raise FiducialConversionError(
            f"output already exists: {names}; pass overwrite=True or --overwrite"
        )

    afni_copy = shutil.which("3dcopy")
    if afni_copy is None:
        raise FiducialConversionError("AFNI command '3dcopy' was not found on PATH")

    try:
        destination.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise FiducialConversionError(
            f"could not create output directory {destination}: {error}"
        ) from error

    if overwrite:
        try:
            for path in existing:
                path.unlink()
        except OSError as error:
            raise FiducialConversionError(
                f"could not remove existing AFNI output {path}: {error}"
            ) from error

    try:
        subprocess.run([afni_copy, str(image), str(prefix)], check=True)
    except (OSError, subprocess.CalledProcessError) as error:
        for path in outputs:
            path.unlink(missing_ok=True)
        raise FiducialConversionError(f"3dcopy failed for {image}: {error}") from error

    generated_brik = brik if brik.is_file() else compressed_brik
    if not head.is_file() or not generated_brik.is_file():
        for path in outputs:
            path.unlink(missing_ok=True)
        raise FiducialConversionError(
            f"3dcopy did not create the expected AFNI dataset at {prefix}+orig"
        )

    try:
        _append_tagset(head, lps_coordinates)
    except Exception:
        for path in (head, brik, compressed_brik):
            if path.exists():
                path.unlink()
        raise

    return generated_brik, head


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="convert_json_fids_to_head",
        description=(
            "Convert a BIDS T1w NIfTI and JSON fiducials to AFNI BRIK/HEAD files."
        ),
    )
    parser.add_argument("t1w_nii", type=Path, help="BIDS T1w image ending in .nii.gz")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="output directory (defaults to the input image directory)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace an existing AFNI output pair",
    )
    args = parser.parse_args()

    try:
        brik, head = convert_json_fids_to_head(
            args.t1w_nii,
            args.output_dir,
            overwrite=args.overwrite,
        )
    except FiducialConversionError as error:
        parser.exit(1, f"{parser.prog}: error: {error}\n")
    print(brik)
    print(head)


if __name__ == "__main__":
    main()
