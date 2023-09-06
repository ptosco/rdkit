""" contains factory class for producing signatures


"""
from __future__ import annotations
import rdkit.Chem.Pharm2D.SigFactory
import typing
from rdkit.DataStructs.cDataStructs import IntSparseIntVect
from rdkit.DataStructs.cDataStructs import LongSparseIntVect
from rdkit.DataStructs.cDataStructs import SparseBitVect
import copy
import numpy
import rdkit.Chem.Pharm2D.Utils
_Shape = typing.Tuple[int, ...]

__all__ = [
    "IntSparseIntVect",
    "LongSparseIntVect",
    "SigFactory",
    "SparseBitVect",
    "Utils",
    "copy",
    "numpy"
]


class SigFactory():
    """
    SigFactory's are used by creating one, setting the relevant
    parameters, then calling the GetSignature() method each time a
    signature is required.
    """
    pass
_verbose = False
