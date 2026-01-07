"""
A standalone bacterial amplicon designing tool (Tuberculosis Optimized Amplicon Sequencing Tool)
"""

# __version__ = '1.5.4'

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("toast-amplicon")
except PackageNotFoundError:
    __version__ = "unknown"

__all__ = ["__version__"]
