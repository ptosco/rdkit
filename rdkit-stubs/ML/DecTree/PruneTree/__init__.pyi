""" Contains functionality for doing tree pruning

"""
from __future__ import annotations
import rdkit.ML.DecTree.PruneTree
import typing
import copy
import numpy
import rdkit.ML.DecTree.CrossValidate
import rdkit.ML.DecTree.DecTree
_Shape = typing.Tuple[int, ...]

__all__ = [
    "CrossValidate",
    "DecTree",
    "MaxCount",
    "PruneTree",
    "copy",
    "numpy"
]


_verbose = 0
