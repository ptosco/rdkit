""" A module for Kier and Hall's EState Descriptors

Unless otherwise noted, all definitions here can be found in:

  L.B. Kier and L.H. Hall _Molecular Structure Description:
  The Electrotopological State"_  Academic Press (1999)

"""
from __future__ import annotations
import rdkit.Chem.EState
import typing
import numpy
import rdkit.Chem
import sys
_Shape = typing.Tuple[int, ...]

__all__ = [
    "AtomTypes",
    "BuildPatts",
    "Chem",
    "EState",
    "EStateIndices",
    "EState_VSA",
    "GetPrincipleQuantumNumber",
    "MaxAbsEStateIndex",
    "MaxEStateIndex",
    "MinAbsEStateIndex",
    "MinEStateIndex",
    "TypeAtoms",
    "esPatterns",
    "numpy",
    "sys"
]


esPatterns = None
