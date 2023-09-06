"""
Supplies a class for working with molecules from databases
"""
from __future__ import annotations
import rdkit.Chem.Suppliers.DbMolSupplier
import typing
from rdkit.Chem.Suppliers.MolSupplier import MolSupplier
import rdkit.Chem
import rdkit.Chem.Suppliers.MolSupplier
import sys

__all__ = [
    "Chem",
    "DbMolSupplier",
    "ForwardDbMolSupplier",
    "MolSupplier",
    "RandomAccessDbMolSupplier",
    "sys",
    "warning"
]


class DbMolSupplier(rdkit.Chem.Suppliers.MolSupplier.MolSupplier):
    """
    new molecules come back with all additional fields from the
    database set in a "_fieldsFromDb" data member
    """
    pass
class ForwardDbMolSupplier(DbMolSupplier, rdkit.Chem.Suppliers.MolSupplier.MolSupplier):
    """
    DbMol supplier supporting only forward iteration


       new molecules come back with all additional fields from the
       database set in a "_fieldsFromDb" data member

     
    """
    pass
class RandomAccessDbMolSupplier(DbMolSupplier, rdkit.Chem.Suppliers.MolSupplier.MolSupplier):
    pass
