from __future__ import annotations
import rdkit.ML.Cluster.murtagh_test
import typing
import numpy
import rdkit.ML.Cluster.Murtagh
_Shape = typing.Tuple[int, ...]

__all__ = [
    "Murtagh",
    "clusters",
    "d",
    "dist",
    "dists",
    "i",
    "j",
    "numpy"
]


clusters: list # value = [<rdkit.ML.Cluster.Clusters.Cluster object>]
d: numpy.ndarray # value = 
"""
array([[10.,  5.],
       [20., 20.],
       [30., 10.],
       [30., 15.],
       [ 5., 10.]])
"""
dist = 650.0
dists: numpy.ndarray # value = array([325., 425., 200., 500., 125.,  25.,  50., 325., 625., 650.])
i = 4
j = 3
