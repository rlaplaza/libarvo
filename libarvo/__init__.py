"""Python library for libarvo."""

from importlib import metadata

from libarvo.arvo import molecular_vs

__all__ = ["molecular_vs"]

# Version
__version__ = metadata.version("libarvo")
