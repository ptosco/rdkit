from __future__ import annotations
import rdkit.Chem.Pharm3D.Pharmacophore
import typing
import numpy
import rdkit.Chem.ChemicalFeatures
import rdkit.Geometry
import rdkit.RDLogger
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ChemicalFeatures",
    "ExplicitPharmacophore",
    "Geometry",
    "Pharmacophore",
    "logger",
    "numpy"
]


class ExplicitPharmacophore():
    """
    this is a pharmacophore with explicit point locations and radii
     
    """
    pass
class Pharmacophore():
    pass
logger: rdkit.RDLogger.logger
