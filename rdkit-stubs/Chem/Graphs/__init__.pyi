""" Python functions for manipulating molecular graphs

In theory much of the functionality in here should be migrating into the
C/C++ codebase.

"""
from __future__ import annotations
import rdkit.Chem.Graphs
import typing
import numpy
import rdkit.Chem
import rdkit.DataStructs
import types
_Shape = typing.Tuple[int, ...]

__all__ = [
    "CharacteristicPolynomial",
    "Chem",
    "DataStructs",
    "numpy",
    "types"
]


