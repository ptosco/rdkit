"""

Here we're implementing the Decision Tree stuff found in Chapter 3 of
Tom Mitchell's Machine Learning Book.

"""
from __future__ import annotations
import rdkit.ML.DecTree
import typing

__all__ = [
    "Tree",
    "TreeUtils"
]


