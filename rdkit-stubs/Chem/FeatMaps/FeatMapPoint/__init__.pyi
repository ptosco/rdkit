from __future__ import annotations
import rdkit.Chem.FeatMaps.FeatMapPoint
import typing
import Boost.Python
import rdkit.Chem.ChemicalFeatures
import rdkit.Chem.rdChemicalFeatures

__all__ = [
    "ChemicalFeatures",
    "FeatMapPoint"
]


class FeatMapPoint(rdkit.Chem.rdChemicalFeatures.FreeChemicalFeature, Boost.Python.instance):
    featDirs = None
    weight = 0.0
    pass
