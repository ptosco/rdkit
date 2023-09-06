from __future__ import annotations
import rdkit.Chem.FastSDMolSupplier
import typing
import Boost.Python
import rdkit.Chem
import rdkit.Chem.rdmolfiles
import sys
import warnings

__all__ = [
    "Chem",
    "FastSDMolSupplier",
    "sys",
    "warnings"
]


class FastSDMolSupplier(rdkit.Chem.rdmolfiles.SDMolSupplier, Boost.Python.instance):
    pass
__warningregistry__ = {'version': 73}
