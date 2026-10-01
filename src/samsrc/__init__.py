"""Python launchers and resources for samsrc."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("samsrc")
except PackageNotFoundError:  # Source-tree imports during development.
    __version__ = "5.1.0"

__all__ = ["__version__"]
