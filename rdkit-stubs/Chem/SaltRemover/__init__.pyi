from __future__ import annotations
import rdkit.Chem.SaltRemover
import typing
from rdkit.Chem.rdmolfiles import SDMolSupplier
from rdkit.Chem.rdmolfiles import SmilesMolSupplier
from contextlib import closing
import os
import rdkit.Chem
import rdkit.RDConfig
import re

__all__ = [
    "Chem",
    "InputFormat",
    "RDConfig",
    "SDMolSupplier",
    "SaltRemover",
    "SmilesMolSupplier",
    "closing",
    "namedtuple",
    "os",
    "re"
]


class InputFormat():
    MOL = 'mol'
    SMARTS = 'smarts'
    SMILES = 'smiles'
    pass
class SaltRemover():
    defnFilename = '/scratch/toscopa1/src/rdkit/Data/Salts.txt'
    pass
