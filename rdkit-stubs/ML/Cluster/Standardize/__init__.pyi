""" contains code for standardization of data matrices for clustering


"""
from __future__ import annotations
import rdkit.ML.Cluster.Standardize
import typing
import rdkit.ML.Data.Stats

__all__ = [
    "Stats",
    "StdDev",
    "methods"
]


methods: list # value = [('None', <function <lambda>>, 'No Standardization'), ('Standard Deviation', <function StdDev>, 'Use the standard deviation')]
