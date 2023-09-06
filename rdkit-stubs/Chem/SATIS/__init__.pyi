"""  Functionality for SATIS typing atoms

"""
from __future__ import annotations
import rdkit.Chem.SATIS
import typing
import itertools
import rdkit.Chem
import rdkit.Chem.rdchem

__all__ = [
    "Chem",
    "SATISTypes",
    "aldehydePatt",
    "amidePatt",
    "carboxylPatt",
    "carboxylatePatt",
    "esterPatt",
    "itertools",
    "ketonePatt",
    "specialCases"
]


aldehydePatt: rdkit.Chem.rdchem.Mol
amidePatt: rdkit.Chem.rdchem.Mol
carboxylPatt: rdkit.Chem.rdchem.Mol
carboxylatePatt: rdkit.Chem.rdchem.Mol
esterPatt: rdkit.Chem.rdchem.Mol
ketonePatt: rdkit.Chem.rdchem.Mol
specialCases: tuple # value = ((<rdkit.Chem.rdchem.Mol object>, 97), (<rdkit.Chem.rdchem.Mol object>, 96), (<rdkit.Chem.rdchem.Mol object>, 98), (<rdkit.Chem.rdchem.Mol object>, 95), (<rdkit.Chem.rdchem.Mol object>, 94), (<rdkit.Chem.rdchem.Mol object>, 93))
