""" ID3 Decision Trees

  contains an implementation of the ID3 decision tree algorithm
  as described in Tom Mitchell's book "Machine Learning"

  It relies upon the _Tree.TreeNode_ data structure (or something
    with the same API) defined locally to represent the trees

"""
from __future__ import annotations
import rdkit.ML.DecTree.ID3
import typing
import numpy
import rdkit.ML.DecTree.DecTree
import rdkit.ML.InfoTheory.entropy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "CalcTotalEntropy",
    "DecTree",
    "GenVarTable",
    "ID3",
    "ID3Boot",
    "entropy",
    "numpy"
]


