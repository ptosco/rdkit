""" code for dealing with forests (collections) of decision trees

**NOTE** This code should be obsolete now that ML.Composite.Composite is up and running.

"""
from __future__ import annotations
import rdkit.ML.DecTree.Forest
import typing
import numpy
import pickle
import rdkit.ML.DecTree.CrossValidate
import rdkit.ML.DecTree.PruneTree
_Shape = typing.Tuple[int, ...]

__all__ = [
    "CrossValidate",
    "Forest",
    "PruneTree",
    "numpy",
    "pickle"
]


class Forest():
    """
    a forest of unique decision trees.

        adding an existing tree just results in its count field being incremented
            and the errors being averaged.

        typical usage:

          1) grow the forest with AddTree until happy with it

          2) call AverageErrors to calculate the average error values

          3) call SortTrees to put things in order by either error or count

      
    """
    pass
