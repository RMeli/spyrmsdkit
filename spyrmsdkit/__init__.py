"""
spyrmsdkit
MDAnalysis Kit (MDAKit) for SPyRMSD
"""

# Add imports here
from .spyrmsdkit import canvas

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
