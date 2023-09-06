"""  EState fingerprinting

"""
from __future__ import annotations
import rdkit.Chem.EState.Fingerprinter
import typing
import numpy
import rdkit.Chem.EState.AtomTypes
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AtomTypes",
    "EStateIndices",
    "FingerprintMol",
    "numpy"
]


