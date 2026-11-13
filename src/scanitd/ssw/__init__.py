"""Vendored SSW (Striped Smith-Waterman) package.

The underlying SSW C library and wrapper conventions are MIT-licensed;
see LICENSE-ssw-py.txt for details.
"""

from .alignment import Alignment, AlignmentMgr

__all__ = ["Alignment", "AlignmentMgr"]
