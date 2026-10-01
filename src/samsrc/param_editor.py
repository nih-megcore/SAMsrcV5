"""SAM parameter-file catalog, parser, serializer, and validation helpers."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import re
from typing import Iterable


PROGRAMS = ("sam_cov", "sam_wts", "sam_3d", "sam_ers")


@dataclass(frozen=True)
class ParameterSpec:
    key: str
    category: str
    kind: str
    programs: frozenset[str]
    help: str
    choices: tuple[str, ...] = ()
    browse: str | None = None
    placeholder: str = ""


def _programs(*names: str) -> frozenset[str]:
    return frozenset(names or PROGRAMS)


SPECS = (
    ParameterSpec("%include", "General", "repeat", _programs(),
                  "Read another parameter file at this point.", placeholder="one filename per line"),
    ParameterSpec("Verbose", "General", "flag", _programs(), "Enable verbose output."),
    ParameterSpec("DataSet", "General", "string", _programs(), "MEG dataset name or path.", browse="directory"),
    ParameterSpec("PDFName", "General", "string", _programs(), "BTi/4D PDF filename.", browse="file"),
    ParameterSpec("InputSAMDirectory", "General", "string", _programs(),
                  "SAM root used to read covariance matrices and weights.", browse="directory"),
    ParameterSpec("OutputSAMDirectory", "General", "string", _programs(),
                  "SAM root used to write products and run logs.", browse="directory"),

    ParameterSpec("Marker", "Markers and Time", "repeat", _programs(),
                  "Marker name, start, end, SUM flag, and optional covariance name.",
                  placeholder="name start end TRUE|FALSE [covariance-name]"),
    ParameterSpec("SegFile", "Markers and Time", "repeat", _programs("sam_cov", "sam_wts", "sam_3d"),
                  "Marker name, segment file, and SUM flag.",
                  placeholder="name filename TRUE|FALSE"),
    ParameterSpec("DataSegment", "Markers and Time", "range", _programs("sam_cov", "sam_ers"),
                  "Processing interval relative to each marker, in seconds.", placeholder="start end"),
    ParameterSpec("Baseline", "Markers and Time", "range", _programs("sam_ers"),
                  "Baseline interval relative to the first marker, in seconds.", placeholder="start end"),
    ParameterSpec("SignSegment", "Markers and Time", "range", _programs("sam_ers"),
                  "Interval used to resolve source polarity, in seconds.", placeholder="start end"),
    ParameterSpec("TimeStep", "Markers and Time", "float", _programs("sam_ers"),
                  "Time between output images, in seconds.", placeholder="seconds"),

    ParameterSpec("XBounds", "Coordinates and MRI", "bounds", _programs("sam_wts"),
                  "Anterior-posterior ROI bounds in centimetres.", placeholder="start end"),
    ParameterSpec("YBounds", "Coordinates and MRI", "bounds", _programs("sam_wts"),
                  "Left-right ROI bounds in centimetres.", placeholder="start end"),
    ParameterSpec("ZBounds", "Coordinates and MRI", "bounds", _programs("sam_wts"),
                  "Inferior-superior ROI bounds in centimetres.", placeholder="start end"),
    ParameterSpec("ImageStep", "Coordinates and MRI", "positive", _programs("sam_wts"),
                  "Image voxel spacing in centimetres.", placeholder="centimetres"),
    ParameterSpec("MRIDirectory", "Coordinates and MRI", "string", _programs("sam_wts"),
                  "Root containing participant MRI files.", browse="directory"),
    ParameterSpec("MRIPattern", "Coordinates and MRI", "string", _programs("sam_wts"),
                  "MRI filename pattern; default is %M/%P/%s.", placeholder="%M/%P/%s"),
    ParameterSpec("PrefixLength", "Coordinates and MRI", "prefix", _programs("sam_wts", "sam_3d", "sam_ers"),
                  "Dataset prefix length or a single delimiter character.", placeholder="integer or delimiter"),
    ParameterSpec("HullName", "Coordinates and MRI", "string", _programs("sam_wts"),
                  "Hull filename; default is hull.shape.", browse="file"),
    ParameterSpec("AtlasName", "Coordinates and MRI", "string", _programs("sam_wts", "sam_3d"),
                  "Atlas containing source coordinates and orientations.", browse="file"),
    ParameterSpec("TargetName", "Coordinates and MRI", "string", _programs("sam_wts"),
                  "File containing discrete target coordinates.", browse="file"),
    ParameterSpec("Extent", "Coordinates and MRI", "positive", _programs("sam_wts"),
                  "Radial extent around targets, in millimetres.", placeholder="millimetres"),
    ParameterSpec("Transform", "Coordinates and MRI", "string", _programs("sam_wts"),
                  "Transform name (without .xfm) read from the SAM input root."),
    ParameterSpec("ImageFormat", "Coordinates and MRI", "imageformat", _programs("sam_wts", "sam_3d", "sam_ers"),
                  "ORIG output, or TLRC followed by voxel resolution in millimetres."),
    ParameterSpec("ImageDirectory", "Coordinates and MRI", "string", _programs("sam_3d", "sam_ers"),
                  "Directory where image files are written.", browse="directory"),

    ParameterSpec("CovBand", "Frequency", "band", _programs(),
                  "Covariance passband in Hz.", placeholder="low high"),
    ParameterSpec("ImageBand", "Frequency", "band", _programs("sam_cov", "sam_3d", "sam_ers"),
                  "Imaging passband in Hz.", placeholder="low high"),
    ParameterSpec("OrientBand", "Frequency", "band", _programs("sam_cov"),
                  "Passband used for orientation covariance in Hz.", placeholder="low high"),
    ParameterSpec("NoiseBand", "Frequency", "band", _programs("sam_cov", "sam_wts", "sam_3d"),
                  "Passband used to estimate sensor noise in Hz.", placeholder="low high"),
    ParameterSpec("SmoothBand", "Frequency", "band", _programs("sam_ers"),
                  "Smoothing passband in Hz.", placeholder="low high"),
    ParameterSpec("FilterType", "Frequency", "choice", _programs("sam_cov", "sam_3d", "sam_ers"),
                  "Filtering implementation.", choices=("FFT", "IIR")),
    ParameterSpec("Notch", "Frequency", "flag", _programs("sam_cov", "sam_3d", "sam_ers"),
                  "Apply mains-frequency notch filters."),
    ParameterSpec("Hz", "Frequency", "positive", _programs("sam_cov", "sam_3d", "sam_ers"),
                  "Electrical mains frequency; default is 60 Hz.", placeholder="50 or 60"),

    ParameterSpec("CovType", "SAM and Imaging", "choice", _programs("sam_3d", "sam_ers"),
                  "Covariance/weight set used by the image analysis.", choices=("GLOBAL", "SUM", "ALL")),
    ParameterSpec("Model", "SAM and Imaging", "model", _programs("sam_wts"),
                  "Forward model: SingleSphere x y z, MultiSphere, or Nolte [order]."),
    ParameterSpec("Order", "SAM and Imaging", "int", _programs("sam_wts"),
                  "Spherical-harmonic order for the Nolte model."),
    ParameterSpec("Mu", "SAM and Imaging", "mu", _programs("sam_wts"),
                  "Regularization as [+]value or *value."),
    ParameterSpec("Pinv", "SAM and Imaging", "int", _programs("sam_cov", "sam_wts"),
                  "Number of smallest dimensions removed by the pseudo-inverse."),
    ParameterSpec("Normalize", "SAM and Imaging", "flag", _programs("sam_wts"),
                  "Normalize SAM weights by projected noise."),
    ParameterSpec("Field", "SAM and Imaging", "flag", _programs("sam_wts"),
                  "Write primary-sensor forward solutions with weights."),
    ParameterSpec("ImageMetric", "SAM and Imaging", "imagemetric", _programs("sam_3d", "sam_ers"),
                  "Output metric. sam_3d supports Power; sam_ers supports Signal or Power."),
    ParameterSpec("Absolute", "SAM and Imaging", "flag", _programs("sam_ers"),
                  "Write absolute ERS voxel values."),
    ParameterSpec("CovName", "SAM and Imaging", "string", _programs("sam_cov", "sam_wts", "sam_3d", "sam_ers"),
                  "Name used when locating covariance products."),
    ParameterSpec("WtsName", "SAM and Imaging", "string", _programs("sam_wts", "sam_3d", "sam_ers"),
                  "Name used when locating weight products."),
    ParameterSpec("OutName", "SAM and Imaging", "string", _programs(),
                  "Name used in output products instead of the parameter filename."),
)

SPEC_BY_KEY = {spec.key: spec for spec in SPECS}
CATEGORIES = tuple(dict.fromkeys(spec.category for spec in SPECS))


def canonical_key(token: str) -> str | None:
    """Resolve parser-style case-insensitive unique abbreviations."""
    if re.fullmatch(r"marker\d+", token, re.IGNORECASE):
        return "Marker"
    lowered = token.casefold()
    matches = [key for key in SPEC_BY_KEY if key.casefold().startswith(lowered)]
    return matches[0] if len(matches) == 1 else None


@dataclass
class DocumentLine:
    raw: str
    key: str | None = None
    arguments: str = ""
    comment: str = ""


@dataclass
class ParameterDocument:
    lines: list[DocumentLine]
    newline: str = "\n"
    source: str | None = None

    @classmethod
    def parse(cls, text: str) -> "ParameterDocument":
        newline = "\r\n" if "\r\n" in text else "\n"
        lines: list[DocumentLine] = []
        for original in text.splitlines():
            stripped = original.lstrip()
            if not stripped or stripped.startswith("#"):
                lines.append(DocumentLine(original))
                continue
            body, separator, comment_text = original.partition("#")
            fields = body.strip().split(None, 1)
            key = canonical_key(fields[0])
            if key is None:
                lines.append(DocumentLine(original))
                continue
            arguments = fields[1].strip() if len(fields) == 2 else ""
            comment = f"#{comment_text}" if separator else ""
            lines.append(DocumentLine(original, key, arguments, comment))
        return cls(lines, newline, text)

    def values(self) -> dict[str, list[str]]:
        values: dict[str, list[str]] = {}
        for line in self.lines:
            if line.key is not None:
                values.setdefault(line.key, []).append(line.arguments)
        return values

    def render(self, updates: dict[str, list[str]]) -> str:
        """Render updates while leaving untouched and unknown lines byte-stable."""
        if not updates and self.source is not None:
            return self.source
        output: list[str] = []
        emitted: set[str] = set()
        for line in self.lines:
            if line.key is None or line.key not in updates:
                output.append(line.raw)
                continue
            if line.key in emitted:
                continue
            emitted.add(line.key)
            replacements = updates[line.key]
            for index, arguments in enumerate(replacements):
                value = line.key if not arguments else f"{line.key} {arguments}"
                if index == 0 and line.comment:
                    value = f"{value}  {line.comment}"
                output.append(value)
        pending = [key for key in SPEC_BY_KEY if key in updates and key not in emitted]
        if pending and output and output[-1] != "":
            output.append("")
        for key in pending:
            for arguments in updates[key]:
                output.append(key if not arguments else f"{key} {arguments}")
        return self.newline.join(output).rstrip() + self.newline


def serialize_new(values: dict[str, list[str]]) -> str:
    return ParameterDocument([]).render(values)


def changed_updates(
    original: dict[str, list[str]], current: dict[str, list[str]]
) -> dict[str, list[str]]:
    keys = set(original) | set(current)
    return {key: current.get(key, []) for key in keys if original.get(key, []) != current.get(key, [])}


def _numbers(arguments: str, count: int) -> list[float]:
    fields = arguments.split()
    if len(fields) != count:
        raise ValueError(f"requires {count} numeric value{'s' if count != 1 else ''}")
    values = [float(field) for field in fields]
    if not all(math.isfinite(value) for value in values):
        raise ValueError("values must be finite")
    return values


def _validate_one(spec: ParameterSpec, arguments: str) -> str | None:
    fields = arguments.split()
    if spec.kind == "flag":
        return None if not fields else "does not take a value"
    if not arguments:
        return "requires a value"
    try:
        if spec.kind in {"float", "positive"}:
            values = _numbers(arguments, 1)
            if spec.kind == "positive" and values[0] <= 0:
                return "must be greater than zero"
        elif spec.kind == "int":
            if len(fields) != 1 or not re.fullmatch(r"[+-]?\d+", fields[0]):
                return "requires one integer"
        elif spec.kind in {"range", "bounds", "band"}:
            values = _numbers(arguments, 2)
            if spec.kind == "bounds" and values[0] > values[1]:
                return "start must not exceed end"
            if spec.kind in {"range", "band"} and values[0] >= values[1]:
                return "start/low value must be less than end/high value"
            if spec.kind == "band" and values[0] < 0:
                return "frequencies must be non-negative"
        elif spec.kind == "choice":
            if len(fields) != 1 or fields[0].upper() not in spec.choices:
                return f"must be one of {', '.join(spec.choices)}"
        elif spec.kind == "prefix":
            if len(fields) != 1:
                return "requires one integer or delimiter character"
            if not re.fullmatch(r"[+-]?\d+", fields[0]) and len(fields[0]) != 1:
                return "must be an integer or one delimiter character"
        elif spec.kind == "mu":
            if len(fields) != 1 or not re.fullmatch(
                r"[+*]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", fields[0]
            ):
                return "must be a number optionally prefixed with + or *"
        elif spec.kind == "imageformat":
            if fields[0].upper() == "ORIG" and len(fields) == 1:
                return None
            if fields[0].upper() == "TLRC" and len(fields) == 2:
                value = _numbers(fields[1], 1)[0]
                return None if value > 0 else "TLRC resolution must be greater than zero"
            return "must be ORIG or TLRC followed by resolution in mm"
        elif spec.kind == "model":
            model = fields[0].casefold()
            if model == "singlesphere" and len(fields) == 4:
                _numbers(" ".join(fields[1:]), 3)
                return None
            if model == "multisphere" and len(fields) == 1:
                return None
            if model == "nolte" and len(fields) in {1, 2}:
                if len(fields) == 2 and not re.fullmatch(r"\d+", fields[1]):
                    return "Nolte order must be an integer"
                return None
            return "must be SingleSphere x y z, MultiSphere, or Nolte [order]"
        elif spec.kind == "imagemetric":
            if fields[0].casefold() not in {"power", "signal"} or len(fields) != 1:
                return "must be Signal or Power"
        elif spec.kind == "repeat":
            if spec.key == "Marker":
                if len(fields) not in {4, 5}:
                    return "requires name, start, end, TRUE|FALSE, and optional covariance name"
                _numbers(" ".join(fields[1:3]), 2)
                if fields[3].upper() not in {"TRUE", "FALSE"}:
                    return "SUM flag must be TRUE or FALSE"
            elif spec.key == "SegFile":
                if len(fields) != 3 or fields[2].upper() not in {"TRUE", "FALSE"}:
                    return "requires marker, filename, and TRUE|FALSE"
            elif len(fields) != 1:
                return "requires one filename"
        elif spec.kind == "string" and len(fields) != 1:
            return "SAM parameter values cannot contain whitespace"
    except ValueError as error:
        return str(error)
    return None


def validate_values(
    values: dict[str, list[str]], profile: str
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    for key, entries in values.items():
        spec = SPEC_BY_KEY.get(key)
        if spec is None:
            continue
        for arguments in entries:
            problem = _validate_one(spec, arguments)
            if problem:
                errors.append(f"{key}: {problem}")

    required = {
        "sam_cov": ("CovBand",),
        "sam_wts": ("CovBand", "Model"),
        "sam_3d": ("CovType", "CovBand", "ImageBand", "ImageMetric"),
        "sam_ers": ("Marker", "CovType", "CovBand", "ImageBand", "SmoothBand", "TimeStep", "ImageMetric"),
    }
    for key in required.get(profile, ()):
        if not values.get(key):
            warnings.append(f"{key} is normally required by {profile}")

    if profile == "sam_wts":
        grid = all(values.get(key) for key in ("XBounds", "YBounds", "ZBounds", "ImageStep"))
        if not grid and not values.get("AtlasName") and not values.get("TargetName"):
            warnings.append("sam_wts normally needs a complete XYZ grid, AtlasName, or TargetName")
        model = values.get("Model", [""])[0].casefold()
        if model.startswith(("nolte", "multisphere")) and not values.get("MRIDirectory"):
            warnings.append("MRIDirectory is required for Nolte and MultiSphere models")
    image_metric = values.get("ImageMetric")
    if (
        profile == "sam_3d"
        and image_metric
        and image_metric[0].casefold() != "power"
    ):
        errors.append("ImageMetric: sam_3d supports only Power")
    if profile == "sam_ers":
        markers = values.get("Marker", [])
        if len(markers) not in {0, 1, 2}:
            errors.append("Marker: sam_ers requires one or two markers")
        covtype = values.get("CovType")
        if covtype and covtype[0].upper() not in {"GLOBAL", "SUM"}:
            errors.append("CovType: sam_ers supports only GLOBAL or SUM")
        metric = values.get("ImageMetric")
        if metric and metric[0].casefold() not in {"signal", "power"}:
            errors.append("ImageMetric: sam_ers supports only Signal or Power")
    cov = values.get("CovBand")
    image = values.get("ImageBand")
    if cov and image:
        try:
            cov_values = _numbers(cov[0], 2)
            image_values = _numbers(image[0], 2)
            if image_values[0] < cov_values[0] or image_values[1] > cov_values[1]:
                errors.append("ImageBand must lie within CovBand")
        except ValueError:
            pass
    return errors, warnings


def managed_parameter_path(root: Path, name: str) -> Path:
    """Return a safe .param target directly beneath root."""
    filename = Path(name).name
    if filename in {"", ".", ".."}:
        raise ValueError("a parameter filename is required")
    if not filename.endswith(".param"):
        filename += ".param"
    return root / filename


def ensure_save_root(home: Path | None = None) -> Path:
    """Create and return the per-user SAM parameter directory."""
    root = (home if home is not None else Path.home()) / "samparams"
    root.mkdir(parents=True, exist_ok=True)
    return root


def lines_to_values(items: Iterable[tuple[str, Iterable[str]]]) -> dict[str, list[str]]:
    return {key: [value.strip() for value in values if value.strip()] for key, values in items}
