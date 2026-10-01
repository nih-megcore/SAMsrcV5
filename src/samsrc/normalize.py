"""Portable implementation of the historical 3dNormalize shell command."""

from __future__ import annotations

import argparse
import gzip
import math
import shutil
import statistics
import subprocess
import sys
from pathlib import Path


def _run(*args: str, capture: bool = False) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(
            args,
            check=True,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
        )
    except FileNotFoundError as error:
        raise SystemExit(
            f"3dNormalize requires AFNI command {args[0]!r} on PATH"
        ) from error
    except subprocess.CalledProcessError as error:
        if capture and error.stderr:
            sys.stderr.write(error.stderr)
        raise SystemExit(error.returncode) from error


def _values(image: str, mask_index: int) -> list[float]:
    result = _run(
        "3dmaskave",
        "-quiet",
        "-mask",
        "SELF",
        "-mindex",
        str(mask_index),
        "-dump",
        image,
        capture=True,
    )
    values = []
    for line in result.stdout.splitlines():
        if line.startswith("+++"):
            continue
        for field in line.split():
            try:
                values.append(float(field))
            except ValueError:
                pass
    if not values:
        raise SystemExit(f"3dNormalize: {image}: no non-zero voxels")
    return values


def _view(input_name: str) -> tuple[str, bool]:
    path = Path(input_name.split("[", 1)[0].split("<", 1)[0])
    name = str(path)
    for suffix in (".gz", ".BRIK", ".HEAD"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    compressed = path.name.endswith(".BRIK.gz") or Path(f"{name}.BRIK.gz").is_file()
    if name.endswith("+acpc"):
        return "acpc", compressed
    if name.endswith("+tlrc"):
        return "tlrc", compressed
    return "orig", compressed


def _tail_sd(values: list[float]) -> float:
    """Match 1dstats on the historical value/-value tail expansion."""
    return math.sqrt(2 * sum(value * value for value in values) / (2 * len(values) - 1))


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="3dNormalize",
        description="Scale an AFNI or SAM volume by its standard deviation.",
    )
    parser.add_argument("-z", action="store_true", help="remove the mean")
    parser.add_argument("-i", action="store_true", help="scale tails independently")
    parser.add_argument("-m", type=int, default=0, metavar="SUBBRIK")
    parser.add_argument("-v", action="store_true")
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()
    if args.z and args.i:
        parser.error("-i and -z cannot be used together")

    values = _values(args.input, args.m)
    view, compressed = _view(args.input)
    output_prefix = args.output
    for candidate in Path(output_prefix).parent.glob(
        f"{Path(output_prefix).name}+{view}*"
    ):
        candidate.unlink()

    if args.i:
        negative = [value for value in values if value < 0]
        positive = [value for value in values if value > 0]
        if not negative or not positive:
            raise SystemExit("3dNormalize: -i requires positive and negative voxels")
        sd_negative = _tail_sd(negative)
        sd_positive = _tail_sd(positive)
        if args.v:
            print(
                f"3dNormalize: {args.input}: -sd is {sd_negative}, "
                f"+sd is {sd_positive}",
                file=sys.stderr,
            )
        expression = (
            f"(isnegative(a)*(a/{sd_negative}))"
            f"+(ispositive(a)*(a/{sd_positive}))"
        )
    else:
        mean = statistics.fmean(values)
        sd = statistics.stdev(values)
        if sd == 0:
            raise SystemExit("3dNormalize: standard deviation is zero")
        if args.v:
            print(
                f"3dNormalize: {args.input}: mean is {mean}, sd is {sd}",
                file=sys.stderr,
            )
        if abs(mean) >= sd:
            relation = "sd < mean" if mean >= sd else "mean < -sd"
            print(f"3dNormalize: {args.input}: {relation}", file=sys.stderr)
        expression = f"bool(a)*(a-({mean}))/{sd}" if args.z else f"bool(a)*(a/{sd})"

    _run("3dcalc", "-a", args.input, "-prefix", output_prefix, "-expr", expression)
    brick = Path(f"{output_prefix}+{view}.BRIK")
    if compressed and brick.exists():
        with brick.open("rb") as source, gzip.open(f"{brick}.gz", "wb") as target:
            shutil.copyfileobj(source, target)
        brick.unlink()


if __name__ == "__main__":
    main()
