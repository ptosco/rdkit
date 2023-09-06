""" code for dealing with resemblance (metric) matrices

    Here's how the matrices are stored:

     '[(0,1),(0,2),(1,2),(0,3),(1,3),(2,3)...]  (row,col), col>row'

     or, alternatively the matrix can be drawn, with indices as:

       || - || 0 || 1 || 3
       || - || - || 2 || 4
       || - || - || - || 5
       || - || - || - || -

     the index of a given (row,col) pair is:
       '(col*(col-1))/2 + row'

"""
from __future__ import annotations
import rdkit.ML.Cluster.Resemblance
import typing
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "CalcMetricMatrix",
    "EuclideanDistance",
    "FindMinValInList",
    "ShowMetricMat",
    "methods",
    "numpy"
]


methods: list # value = [('Euclidean', <function EuclideanDistance>, 'Euclidean Distance')]
