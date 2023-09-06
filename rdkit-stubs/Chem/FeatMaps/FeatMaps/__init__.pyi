from __future__ import annotations
import rdkit.Chem.FeatMaps.FeatMaps
import typing
from rdkit.Chem.FeatMaps.FeatMapPoint import FeatMapPoint
import math

__all__ = [
    "FeatDirScoreMode",
    "FeatMap",
    "FeatMapParams",
    "FeatMapPoint",
    "FeatMapScoreMode",
    "math"
]


class FeatDirScoreMode():
    DotFullRange = 1
    DotPosRange = 2
    Ignore = 0
    pass
class FeatMap():
    dirScoreMode = 0
    params = {}
    scoreMode = 0
    pass
class FeatMapParams():
    """
    one of these should be instantiated for each
     feature type in the feature map
     
    """
    class FeatProfile():
        """
        scoring profile of the feature 
        """
        Box = 2
        Gaussian = 0
        Triangle = 1
        pass
    featProfile = 0
    radius = 2.5
    width = 1.0
    pass
class FeatMapScoreMode():
    All = 0
    Best = 2
    Closest = 1
    pass
