""" Atom-based calculation of LogP and MR using Crippen's approach


    Reference:
      S. A. Wildman and G. M. Crippen *JCICS* _39_ 868-873 (1999)


"""
from __future__ import annotations
import rdkit.Chem.Crippen
import typing
import numpy
import os
import rdkit.Chem
import rdkit.Chem.rdMolDescriptors
import rdkit.RDConfig
_Shape = typing.Tuple[int, ...]

__all__ = [
    "Chem",
    "MolLogP",
    "MolMR",
    "RDConfig",
    "defaultPatternFileName",
    "numpy",
    "os",
    "rdMolDescriptors"
]


_patternOrder = []
_smartsPatterns = {}
defaultPatternFileName = '/scratch/toscopa1/src/rdkit/Data/Crippen.txt'
MolLogP = rdkit.Chem.Crippen.<lambda>
MolMR = rdkit.Chem.Crippen.<lambda>
