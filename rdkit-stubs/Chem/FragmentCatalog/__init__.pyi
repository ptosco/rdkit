from __future__ import annotations
import rdkit.Chem.FragmentCatalog
import typing
from rdkit.Chem.rdfragcatalog import FragCatGenerator
from rdkit.Chem.rdfragcatalog import FragCatParams
from rdkit.Chem.rdfragcatalog import FragCatalog
from rdkit.Chem.rdfragcatalog import FragFPGenerator
import rdkit.Chem
import sys

__all__ = [
    "BitGainsInfo",
    "BuildAdjacencyList",
    "Chem",
    "FragCatGenerator",
    "FragCatParams",
    "FragCatalog",
    "FragFPGenerator",
    "GetMolsMatchingBit",
    "ProcessGainsFile",
    "message",
    "sys"
]


class BitGainsInfo():
    description = ''
    gain = 0.0
    id = -1
    nPerClass = None
    pass
