""" Functionality for ranking bits using info gains

 **Definitions used in this module**

    - *sequence*: an object capable of containing other objects which supports
      __getitem__() and __len__().  Examples of these include lists, tuples, and
      Numeric arrays.

    - *IntVector*: an object containing integers which supports __getitem__() and
       __len__(). Examples include lists, tuples, Numeric Arrays, and BitVects.


 **NOTE**: Neither *sequences* nor *IntVectors* need to support item assignment.
   It is perfectly acceptable for them to be read-only, so long as they are
   random-access.

"""
from __future__ import annotations
import rdkit.ML.InfoTheory.BitRank
import typing
import numpy
import rdkit.ML.InfoTheory.entropy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AnalyzeSparseVects",
    "CalcInfoGains",
    "FormCounts",
    "RankBits",
    "SparseRankBits",
    "entropy",
    "numpy"
]


