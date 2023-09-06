""" Implementation of the clustering algorithm published in:
  Butina JCICS 39 747-750 (1999)

"""
from __future__ import annotations
import rdkit.ML.Cluster.Butina
import typing
import numpy
import rdkit.RDLogger
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ClusterData",
    "EuclideanDist",
    "RDLogger",
    "logger",
    "numpy"
]


logger: rdkit.RDLogger.logger
