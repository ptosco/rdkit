from __future__ import annotations
import rdkit.Chem.FeatMaps.FeatMapParser
import typing
import rdkit.Chem.FeatMaps.FeatMapPoint
import rdkit.Chem.FeatMaps.FeatMaps
import rdkit.Geometry
import re

__all__ = [
    "FeatMapParseError",
    "FeatMapParser",
    "FeatMapPoint",
    "FeatMaps",
    "Geometry",
    "re"
]


class FeatMapParseError(ValueError, Exception, BaseException):
    pass
class FeatMapParser():
    data = None
    pass
