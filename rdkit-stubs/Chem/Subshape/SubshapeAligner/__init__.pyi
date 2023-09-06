from __future__ import annotations
import rdkit.Chem.Subshape.SubshapeAligner
import typing
import numpy
import rdkit.Chem
import rdkit.Chem.Subshape.SubshapeObjects
import rdkit.Geometry
import rdkit.Numerics.rdAlignment
import rdkit.RDLogger
_Shape = typing.Tuple[int, ...]

__all__ = [
    "Alignment",
    "Chem",
    "ClusterAlignments",
    "Geometry",
    "GetShapeShapeDistance",
    "RDLogger",
    "SubshapeAligner",
    "SubshapeAlignment",
    "SubshapeDistanceMetric",
    "SubshapeObjects",
    "TransformMol",
    "logger",
    "numpy"
]


class SubshapeAligner():
    coarseGridToleranceMult = 1.0
    dirThresh = 2.6
    distMetric = 1
    edgeTol = 6.0
    medGridToleranceMult = 1.0
    numFeatThresh = 3
    shapeDistTol = 0.2
    triangleRMSTol = 1.0
    pass
class SubshapeAlignment():
    alignedConfId = -1
    dirMatch = 0.0
    queryTri = None
    shapeDist = 0.0
    targetTri = None
    transform = None
    triangleSSD = None
    pass
class SubshapeDistanceMetric():
    PROTRUDE = 1
    TANIMOTO = 0
    pass
logger: rdkit.RDLogger.logger
