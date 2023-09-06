from __future__ import annotations
import rdkit.Chem.Subshape.SubshapeObjects
import typing

__all__ = [
    "DisplaySubshape",
    "DisplaySubshapeSkeleton",
    "ShapeWithSkeleton",
    "SkeletonPoint",
    "SubshapeShape"
]


class ShapeWithSkeleton():
    grid = None
    skelPts = None
    pass
class SkeletonPoint():
    featmapFeatures = None
    fracVol = 0.0
    location = None
    molFeatures = None
    shapeDirs = None
    shapeMoments = None
    pass
class SubshapeShape():
    featMap = None
    keyFeat = None
    shapes = None
    pass
