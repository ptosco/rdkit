""" Supplies an abstract class for working with sequences of molecules

"""
from __future__ import annotations
import rdkit.Chem.Suppliers.MolSupplier
import typing

__all__ = [
    "MolSupplier"
]


class MolSupplier():
    """
    we must, at minimum, support forward iteration

     
    """
    pass
