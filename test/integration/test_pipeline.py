import json
import math
import os
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "test" / "fixtures"
EXPECTED = json.loads((FIXTURES / "expected.json").read_text())
DATA_ROOT = Path(os.environ.get("TEST_DATA_DIR", ROOT / ".test-data"))
RESULTS_ROOT = Path(os.environ.get("TEST_RESULTS_DIR", ROOT / ".test-results"))


def run(command, cwd, name):
    result = subprocess.run(
        [str(value) for value in command],
        cwd=cwd,
        env={**os.environ, "PWD": str(cwd)},
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    (RESULTS_ROOT / f"{name}.log").write_text(result.stdout)
    assert result.returncode == 0, (
        f"command failed ({result.returncode}): {' '.join(map(str, command))}\n"
        f"{result.stdout}"
    )
    return result.stdout


def afni_info(dataset, *options):
    output = run(["3dinfo", *options, dataset], Path(dataset).parent, "afni-info")
    return output.strip().replace("|", " ").split()


def brick_stats(dataset):
    output = run(
        ["3dBrickStat", "-non-zero", "-mean", "-stdev", dataset],
        Path(dataset).parent,
        f"brick-stat-{Path(dataset).stem}",
    )
    values = [float(value) for value in output.split()[-2:]]
    assert all(math.isfinite(value) for value in values)
    return values


def copy_source_data(work):
    mri_source = DATA_ROOT / "MRI"
    dataset_source = (
        DATA_ROOT / "20010101" / "ABABABAB_airpuff_20010101_001.ds"
    )
    assert mri_source.is_dir(), f"fixture directory missing: {mri_source}"
    assert dataset_source.is_dir(), f"fixture directory missing: {dataset_source}"

    mri_work = work / "MRI"
    mri_work.mkdir(parents=True)
    for suffix in ("HEAD", "BRIK.gz"):
        shutil.copy2(mri_source / f"ABABABAB_refaced+orig.{suffix}", mri_work)

    data_work = work / "20010101"
    data_work.mkdir(parents=True)
    dataset_work = data_work / dataset_source.name
    shutil.copytree(dataset_source, dataset_work)
    return mri_work, data_work, dataset_work


def validate_mri(mri_work):
    dataset = mri_work / "ABABABAB_refaced+orig.HEAD"
    expected = EXPECTED["mri"]

    dimensions = [int(value) for value in afni_info(dataset, "-n4")]
    assert dimensions == expected["dimensions"]
    deltas = [float(value) for value in afni_info(dataset, "-di", "-dj", "-dk")]
    assert deltas == pytest.approx(expected["delta"], abs=1e-5)
    origin = [float(value) for value in afni_info(dataset, "-oi", "-oj", "-ok")]
    assert origin == pytest.approx(expected["origin"], abs=1e-4)
    assert afni_info(dataset, "-orient")[0] == expected["orientation"]
    assert afni_info(dataset, "-space")[0] == expected["space"]
    assert afni_info(dataset, "-datum")[0] == expected["datum"]
    minimum = float(run(["3dBrickStat", "-min", dataset], mri_work, "mri-min").split()[-1])
    maximum = float(run(["3dBrickStat", "-max", dataset], mri_work, "mri-max").split()[-1])
    assert minimum == pytest.approx(expected["minimum"])
    assert maximum == pytest.approx(expected["maximum"])
    head_text = dataset.read_text(errors="replace")
    assert "TAGSET_NUM" in head_text
    for label in expected["tag_labels"]:
        assert label in head_text

    normalized = mri_work / "normalized"
    run(["3dNormalize", "-z", dataset, normalized], mri_work, "afni-normalize")
    normalized_head = mri_work / "normalized+orig.HEAD"
    assert normalized_head.exists()
    mean, standard_deviation = brick_stats(normalized_head)
    assert mean == pytest.approx(0.0, abs=1e-4)
    assert standard_deviation == pytest.approx(1.0, rel=1e-3)
    return dataset


def validate_hull(mri_work, dataset):
    run(["orthohull", "-i0", dataset], mri_work, "orthohull")
    expected_files = (
        "ortho+orig.HEAD",
        "ortho+orig.BRIK",
        "mask+orig.HEAD",
        "mask+orig.BRIK",
        "ortho.ply",
        "ortho_brainhull.ply",
        "ortho_innerskull.ply",
        "ortho_outerskull.ply",
        "SpharmDeco.ply",
        "hull.shape",
        "multisphere.shape",
        "multisphere.shape_info",
    )
    for relative in expected_files:
        path = mri_work / relative
        assert path.is_file() and path.stat().st_size > 0, f"missing output: {path}"

    lines = (mri_work / "hull.shape").read_text().splitlines()
    vertex_count = int(lines[0])
    assert vertex_count == EXPECTED["hull"]["vertices"]
    vertices = np.loadtxt(lines[1 : vertex_count + 1])
    assert vertices.shape == (vertex_count, 6)
    assert np.isfinite(vertices).all()
    assert np.max(np.abs(vertices[:, :3])) < 0.2
    normal_lengths = np.linalg.norm(vertices[:, 3:], axis=1)
    assert normal_lengths == pytest.approx(np.ones(vertex_count), abs=2e-3)


def hull_vertices(path):
    lines = path.read_text().splitlines()
    vertex_count = int(lines[0])
    vertices = np.loadtxt(lines[1 : vertex_count + 1])
    assert vertices.shape == (vertex_count, 6)
    assert np.isfinite(vertices).all()
    return vertices


def validate_ctf_metadata(dataset_work):
    output = run(
        [ROOT / "test" / "bin" / "inspect_ctf", dataset_work],
        dataset_work.parent,
        "inspect-ctf",
    )
    actual = json.loads(output)
    expected = EXPECTED["ctf"]
    for key in ("epochs", "primary_channels", "samples", "markers"):
        assert actual[key] == expected[key]
    assert actual["sample_rate"] == pytest.approx(expected["sample_rate"])
    assert {
        "stim": actual["stim"],
        "missingstim": actual["missingstim"],
    } == expected["marker_counts"]


def read_covariance(path):
    header_format = "<i256s256si4d4i"
    header_size = struct.calcsize(header_format)
    data = path.read_bytes()
    identity = data[:8]
    header = struct.unpack_from(header_format, data, 8)
    channels = header[3]
    channel_index_offset = 8 + header_size
    matrix_offset = channel_index_offset + channels * 4
    assert len(data) == matrix_offset + channels * channels * 8
    covariance = np.frombuffer(
        data,
        dtype="<f8",
        offset=matrix_offset,
        count=channels * channels,
    ).reshape(channels, channels)
    return {
        "identity": identity,
        "version": header[0],
        "set_name": header[1].split(b"\0", 1)[0].decode(),
        "hp": header[4],
        "lp": header[5],
        "type": header[10],
        "segments": header[8],
        "channels": channels,
        "samples": header[9],
        "covariance": covariance,
    }


def validate_covariances(data_work, dataset_work):
    run(
        ["sam_cov", "-r", dataset_work.name, "-m", FIXTURES / "airpuff.param", "-v"],
        data_work,
        "sam-cov",
    )
    covariance_dir = dataset_work / "SAM" / "airpuff,5-70Hz"
    assert covariance_dir.is_dir()
    for name, expected in EXPECTED["covariances"].items():
        path = covariance_dir / f"{name}.cov"
        assert path.is_file()
        actual = read_covariance(path)
        assert actual["identity"] == b"SAMCOVAR"
        assert actual["version"] > 0
        assert actual["set_name"] == dataset_work.stem
        assert actual["hp"] == pytest.approx(5.0)
        assert actual["lp"] == pytest.approx(70.0)
        assert actual["channels"] == EXPECTED["ctf"]["primary_channels"]
        assert actual["type"] == expected["type"]
        assert actual["segments"] == expected["segments"]
        assert actual["samples"] == expected["samples"]
        covariance = actual["covariance"]
        assert np.isfinite(covariance).all()
        assert covariance == pytest.approx(covariance.T, rel=1e-10, abs=1e-12)
        assert np.all(np.diag(covariance) > 0)


def validate_weights_and_images(work, mri_work, data_work, dataset_work):
    mri_root = work / "subjects"
    subject_dir = mri_root / "ABABABAB"
    subject_dir.mkdir(parents=True)
    shutil.copy2(mri_work / "hull.shape", subject_dir)

    weights_parameter = work / "weights.param"
    weights_parameter.write_text(
        (FIXTURES / "weights.param.in").read_text().replace(
            "@@MRI_DIRECTORY@@", str(mri_root)
        )
    )
    run(
        [
            "sam_wts",
            "-r",
            dataset_work.name,
            "-m",
            weights_parameter,
            "-C",
            "airpuff",
            "-v",
        ],
        data_work,
        "sam-wts",
    )

    weights_dir = dataset_work / "SAM" / "airpuff,5-70Hz"
    weights = weights_dir / "Global.nii"
    condition = weights_dir / "GlobalCN.dat"
    assert weights.is_file()
    assert [int(value) for value in afni_info(weights, "-n4")] == [4, 4, 4, 272]
    condition_values = np.loadtxt(condition)
    assert condition_values.shape == (64,)
    assert np.isfinite(condition_values).all()
    assert (weights_dir / "Global_Noise").is_file()

    run(
        [
            "sam_3d",
            "-r",
            dataset_work.name,
            "-m",
            FIXTURES / "image.param",
            "-W",
            "airpuff",
            "-N",
            "integration",
            "-v",
        ],
        data_work,
        "sam-3d",
    )
    image_dir = data_work / "images"
    summary = {}
    for statistic, expected in EXPECTED["images"].items():
        path = image_dir / f"ABABABAB,integration,stim,3D_PWR,{statistic}.nii"
        assert path.is_file()
        assert [int(value) for value in afni_info(path, "-n4")] == [4, 4, 4, 1]
        deltas = [float(value) for value in afni_info(path, "-di", "-dj", "-dk")]
        assert deltas == pytest.approx([10.0, 10.0, -10.0], abs=1e-5)
        assert afni_info(path, "-orient")[0] == "IRP"
        mean, standard_deviation = brick_stats(path)
        assert mean == pytest.approx(expected["mean"], rel=0.05)
        assert standard_deviation == pytest.approx(expected["stdev"], rel=0.05)
        summary[statistic] = {"mean": mean, "stdev": standard_deviation}
        shutil.copy2(path, RESULTS_ROOT / path.name)

    shutil.copy2(mri_work / "hull.shape", RESULTS_ROOT / "hull.shape")
    (RESULTS_ROOT / "pipeline-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )


def validate_decoupled_sam_directories(work, mri_work, data_work, dataset_work):
    covariance_root = work / "external" / "nested" / "covariances"
    weights_root = work / "external" / "weights"
    images_root = work / "external" / "images"

    run(
        [
            "sam_cov",
            "-r",
            dataset_work.name,
            "--Marker",
            "stim",
            "-0.10",
            "0.30",
            "TRUE",
            "--CovBand",
            "5",
            "70",
            "--FilterType",
            "FFT",
            "--Notch",
            "-o_SAMdir",
            covariance_root,
        ],
        data_work,
        "sam-cov-decoupled",
    )
    covariance_dir = covariance_root / "cmdline,5-70Hz"
    assert (covariance_dir / "Global.cov").is_file()
    assert (covariance_root / "sam_cov.param").is_file()
    assert not (dataset_work / "SAM").exists()

    mri_root = work / "decoupled-subjects"
    subject_dir = mri_root / "ABABABAB"
    subject_dir.mkdir(parents=True)
    shutil.copy2(mri_work / "hull.shape", subject_dir)
    run(
        [
            "sam_wts",
            "-r",
            dataset_work.name,
            "--CovBand",
            "5",
            "70",
            "--XBounds",
            "-1",
            "1",
            "--YBounds",
            "-1",
            "1",
            "--ZBounds",
            "5",
            "7",
            "--ImageStep",
            "1",
            "--MRIDirectory",
            mri_root,
            "--Model",
            "Nolte",
            "--Order",
            "8",
            "--ImageFormat",
            "ORIG",
            "-i_SAMdir",
            covariance_root,
            "-o_SAMdir",
            weights_root,
        ],
        data_work,
        "sam-wts-decoupled",
    )
    weights_dir = weights_root / "cmdline,5-70Hz"
    assert (weights_dir / "Global.nii").is_file()
    assert (weights_dir / "GlobalCN.dat").is_file()
    assert (weights_dir / "Global_Noise").is_file()
    assert (weights_root / "sam_wts.param").is_file()
    assert not (covariance_dir / "Global.nii").exists()
    assert not (dataset_work / "SAM").exists()

    epi_parameter = work / "epi.param"
    epi_parameter.write_text("CovBand 5 70\nImageBand 5 70\nTimeInt 2\n")
    epi_root = work / "external" / "epi"
    run(
        [
            "sam_epi",
            "-r",
            dataset_work.name,
            "-m",
            epi_parameter,
            "-W",
            "cmdline",
            "-N",
            "integration",
            "-i_SAMdir",
            weights_root,
            "-o_SAMdir",
            epi_root,
        ],
        data_work,
        "sam-epi-decoupled",
    )
    epi_image = epi_root / "Image" / f"{dataset_work.stem},integration.nii"
    assert epi_image.is_file()
    assert [int(value) for value in afni_info(epi_image, "-n4")] == [4, 4, 4, 1]
    epi_bytes = epi_image.read_bytes()
    epi_offset = int(struct.unpack_from("<f", epi_bytes, 108)[0])
    epi_values = np.frombuffer(epi_bytes, dtype="<f4", count=64, offset=epi_offset)
    assert epi_values.shape == (64,)
    assert np.isfinite(epi_values).all()

    ers_root = work / "external" / "ers"
    run(
        [
            "sam_ers",
            "-r",
            dataset_work.name,
            "--Marker",
            "stim",
            "-0.10",
            "0.30",
            "TRUE",
            "--CovBand",
            "5",
            "70",
            "--ImageBand",
            "5",
            "70",
            "--SmoothBand",
            "0",
            "20",
            "--FilterType",
            "FFT",
            "--CovType",
            "GLOBAL",
            "--ImageMetric",
            "Signal",
            "--TimeStep",
            "0.02",
            "-i_SAMdir",
            weights_root,
            "-o_SAMdir",
            ers_root,
        ],
        data_work,
        "sam-ers-decoupled",
    )
    assert (ers_root / "ABABABAB,cmdline,stim,MOM,ERS.nii").is_file()
    assert (ers_root / "sam_ers.param").is_file()
    assert not (dataset_work / "SAM").exists()

    run(
        [
            "sam_3d",
            "-r",
            dataset_work.name,
            "--Marker",
            "stim",
            "-0.10",
            "0.30",
            "TRUE",
            "--CovBand",
            "5",
            "70",
            "--ImageBand",
            "5",
            "70",
            "--FilterType",
            "FFT",
            "--Notch",
            "--CovType",
            "GLOBAL",
            "--ImageMetric",
            "Power",
            "-i_SAMdir",
            weights_root,
            "-o_SAMdir",
            images_root,
        ],
        data_work,
        "sam-3d-decoupled",
    )
    image = images_root / "ABABABAB,cmdline,stim,3D_PWR,Mean.nii"
    assert image.is_file()
    assert (images_root / "sam_3d.param").is_file()
    assert not (dataset_work / "SAM").exists()

    explicit_images = work / "explicit-images"
    explicit_sam_root = work / "external" / "explicit-sam-output"
    explicit_parameter = work / "explicit-image.param"
    explicit_parameter.write_text(
        (FIXTURES / "image.param").read_text().replace(
            "ImageDirectory images", f"ImageDirectory {explicit_images}"
        )
        + f"InputSAMDirectory {weights_root}\n"
        + f"OutputSAMDirectory {explicit_sam_root}\n"
    )
    run(
        [
            "sam_3d",
            "-r",
            dataset_work.name,
            "-m",
            explicit_parameter,
            "-W",
            "cmdline",
            "-N",
            "explicit",
        ],
        data_work,
        "sam-3d-explicit-images",
    )
    assert (explicit_images / "ABABABAB,explicit,stim,3D_PWR,Mean.nii").is_file()
    assert not list(explicit_sam_root.glob("*.nii"))
    assert (explicit_sam_root / "sam_3d.param").is_file()


@pytest.mark.integration
def test_afni_ctf_pipeline(tmp_path):
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    mri_work, data_work, dataset_work = copy_source_data(tmp_path)
    dataset = validate_mri(mri_work)
    validate_hull(mri_work, dataset)
    validate_ctf_metadata(dataset_work)
    validate_decoupled_sam_directories(tmp_path, mri_work, data_work, dataset_work)
    validate_covariances(data_work, dataset_work)
    validate_weights_and_images(tmp_path, mri_work, data_work, dataset_work)


@pytest.mark.slow
def test_full_brain_5mm_beamformer(tmp_path):
    RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
    slow_results = RESULTS_ROOT / "slow"
    slow_results.mkdir(parents=True, exist_ok=True)

    mri_work, data_work, dataset_work = copy_source_data(tmp_path)
    dataset = mri_work / "ABABABAB_refaced+orig.HEAD"
    validate_hull(mri_work, dataset)
    validate_ctf_metadata(dataset_work)
    validate_covariances(data_work, dataset_work)

    mri_root = tmp_path / "subjects"
    subject_dir = mri_root / "ABABABAB"
    subject_dir.mkdir(parents=True)
    hull = mri_work / "hull.shape"
    shutil.copy2(hull, subject_dir)

    vertices = hull_vertices(hull)
    minimum = vertices[:, :3].min(axis=0) * 100.0
    maximum = vertices[:, :3].max(axis=0) * 100.0
    bounds = list(zip(minimum, maximum))

    weights_parameter = tmp_path / "full5mm.param"
    weights_parameter.write_text(
        "CovBand 5 70\n"
        f"XBounds {bounds[0][0]:.8f} {bounds[0][1]:.8f}\n"
        f"YBounds {bounds[1][0]:.8f} {bounds[1][1]:.8f}\n"
        f"ZBounds {bounds[2][0]:.8f} {bounds[2][1]:.8f}\n"
        "ImageStep 0.5\n"
        f"MRIDirectory {mri_root}\n"
        "Model Nolte\n"
        "Order 8\n"
    )
    run(
        [
            "sam_wts",
            "-r",
            dataset_work.name,
            "-m",
            weights_parameter,
            "-C",
            "airpuff",
            "-W",
            "full5mm",
            "-v",
        ],
        data_work,
        "slow-sam-wts",
    )

    weights_dir = dataset_work / "SAM" / "full5mm,5-70Hz"
    weights = weights_dir / "Global.nii"
    condition = weights_dir / "GlobalCN.dat"
    assert weights.is_file() and condition.is_file()
    dimensions = [int(value) for value in afni_info(weights, "-n4")]
    assert all(size > 4 for size in dimensions[:3])
    assert dimensions[3] == EXPECTED["ctf"]["primary_channels"]
    deltas = [float(value) for value in afni_info(weights, "-di", "-dj", "-dk")]
    assert np.abs(deltas) == pytest.approx([5.0, 5.0, 5.0], abs=1e-5)
    condition_values = np.loadtxt(condition)
    assert condition_values.size == math.prod(dimensions[:3])
    assert np.isfinite(condition_values).all()

    covariance_noise = dataset_work / "SAM" / "airpuff,5-70Hz" / "Global_Noise"
    assert covariance_noise.is_file()
    shutil.copy2(covariance_noise, weights_dir / "Global_Noise")

    image_parameter = tmp_path / "full5mm-image.param"
    image_parameter.write_text(
        (FIXTURES / "image.param").read_text().replace(
            "ImageDirectory images", f"ImageDirectory {slow_results}"
        )
    )
    run(
        [
            "sam_3d",
            "-r",
            dataset_work.name,
            "-m",
            image_parameter,
            "-W",
            "full5mm",
            "-N",
            "full5mm",
            "-v",
        ],
        data_work,
        "slow-sam-3d",
    )

    summary = {
        "bounds_cm": {
            axis: [float(low), float(high)]
            for axis, (low, high) in zip("xyz", bounds)
        },
        "dimensions": dimensions[:3],
        "voxel_size_mm": [abs(value) for value in deltas],
        "images": {},
    }
    for statistic in ("Mean", "Variance"):
        path = slow_results / f"ABABABAB,full5mm,stim,3D_PWR,{statistic}.nii"
        assert path.is_file()
        assert [int(value) for value in afni_info(path, "-n4")] == dimensions[:3] + [1]
        mean, standard_deviation = brick_stats(path)
        assert mean > 0.0
        assert standard_deviation > 0.0
        summary["images"][statistic] = {
            "mean": mean,
            "stdev": standard_deviation,
        }

    shutil.copy2(hull, slow_results / "hull.shape")
    (slow_results / "pipeline-summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
