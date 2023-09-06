""" Descriptors derived from a molecule's 3D structure

"""
from __future__ import annotations
import rdkit.Chem.Descriptors3D
import typing
import rdkit.Chem.rdMolDescriptors

__all__ = [
    "Asphericity",
    "Eccentricity",
    "InertialShapeFactor",
    "NPR1",
    "NPR2",
    "PMI1",
    "PMI2",
    "PMI3",
    "RadiusOfGyration",
    "SpherocityIndex",
    "rdMolDescriptors"
]


Asphericity = rdkit.Chem.Descriptors3D.<lambda>
Eccentricity = rdkit.Chem.Descriptors3D.<lambda>
InertialShapeFactor = rdkit.Chem.Descriptors3D.<lambda>
NPR1 = rdkit.Chem.Descriptors3D.<lambda>
NPR2 = rdkit.Chem.Descriptors3D.<lambda>
PMI1 = rdkit.Chem.Descriptors3D.<lambda>
PMI2 = rdkit.Chem.Descriptors3D.<lambda>
PMI3 = rdkit.Chem.Descriptors3D.<lambda>
RadiusOfGyration = rdkit.Chem.Descriptors3D.<lambda>
SpherocityIndex = rdkit.Chem.Descriptors3D.<lambda>
