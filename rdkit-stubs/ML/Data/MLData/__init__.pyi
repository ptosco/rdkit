""" classes to be used to help work with data sets

"""
from __future__ import annotations
import rdkit.ML.Data.MLData
import typing
import copy
import math
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "MLDataSet",
    "MLQuantDataSet",
    "copy",
    "math",
    "numericTypes",
    "numpy"
]


class MLDataSet():
    """
    A data set for holding general data (floats, ints, and strings)

        **Note**
          this is intended to be a read-only data structure
          (i.e. after calling the constructor you cannot touch it)
       
    """
    pass
class MLQuantDataSet(MLDataSet):
    """
    a data set for holding quantized data


         **Note**

           this is intended to be a read-only data structure
           (i.e. after calling the constructor you cannot touch it)

         **Big differences to MLDataSet**

           1) data are stored in a numpy array since they are homogenous

           2) results are assumed to be quantized (i.e. no qBounds entry is required)

       
    """
    pass
numericTypes: tuple # value = (<class 'int'>, <class 'float'>)
