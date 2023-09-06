""" functionality for finding pharmacophore matches in molecules


  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
from __future__ import annotations
import rdkit.Chem.Pharm2D.Matcher
import typing
import rdkit.Chem
import rdkit.Chem.Pharm2D.Utils

__all__ = [
    "Chem",
    "GetAtomsMatchingBit",
    "MatchError",
    "Utils"
]


class MatchError(Exception, BaseException):
    pass
_verbose = 0
