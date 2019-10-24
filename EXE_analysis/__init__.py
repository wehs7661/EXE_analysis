"""
EXE_analysis
EXE_analysis is a Python package of data analysis tools for expanded ensemble (EXE) simulations.
"""

# Add imports here
from .EXE_analysis import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
