""" handles doing cross validation with decision trees

This is, perhaps, a little misleading.  For the purposes of this module,
cross validation == evaluating the accuracy of a tree.


"""
from __future__ import annotations
import rdkit.ML.DecTree.CrossValidate
import typing
import numpy
import rdkit.ML.Data.SplitData
import rdkit.ML.DecTree.ID3
import rdkit.ML.DecTree.randomtest
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ChooseOptimalRoot",
    "CrossValidate",
    "CrossValidationDriver",
    "ID3",
    "SplitData",
    "TestRun",
    "numpy",
    "randomtest"
]


