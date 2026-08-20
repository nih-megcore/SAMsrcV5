"""Python launchers and resources for SAMsrcV5."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("samsrcv5")
except PackageNotFoundError:  # Source-tree imports during development.
    __version__ = "5.1.0"

__all__ = ["__version__"]
