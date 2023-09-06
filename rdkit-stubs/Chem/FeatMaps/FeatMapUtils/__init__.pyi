from __future__ import annotations
import rdkit.Chem.FeatMaps.FeatMapUtils
import typing
import copy
import rdkit.Chem.FeatMaps.FeatMaps

__all__ = [
    "CombineFeatMaps",
    "DirMergeMode",
    "FeatMaps",
    "GetFeatFeatDistMatrix",
    "MergeFeatPoints",
    "MergeMethod",
    "MergeMetric",
    "copy",
    "familiesMatch",
    "feq"
]


class DirMergeMode():
    NoMerge = 0
    Sum = 1
    pass
class MergeMethod():
    Average = 1
    UseLarger = 2
    WeightedAverage = 0
    pass
class MergeMetric():
    Distance = 1
    NoMerge = 0
    Overlap = 2
    pass
