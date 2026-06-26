"""
StarRotator: Simulating a rotation-broadened stellar spectrum during an exoplanet transit event.

This package provides functionality for loading, processing,
and modeling stellar rotation data.

Modules:
--------
- core: Contains the StarRotator class and all high level functions..
- lib: Utility functions for calculations, operations, etc.
"""

from . import lib # Expose important submodules and functions
from .core import StarRotator

# Optional: Define what gets imported with `from starrotator import *`
# __all__ = ["lib", "StarRotator"]