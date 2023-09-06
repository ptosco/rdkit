""" Utilities for data manipulation

**FILE FORMATS:**

 - *.qdat files* contain quantized data suitable for
  feeding to learning algorithms.

  The .qdat file, written by _DecTreeGui_, is structured as follows:

   1) Any number of lines which are ignored.

   2) A line containing the string 'Variable Table'

      any number of variable definitions in the format:

      '# Variable_name [quant_bounds]'

        where '[quant_bounds]' is a list of the boundaries used for quantizing
         that variable.  If the variable is inherently integral (i.e. not
         quantized), this can be an empty list.

   3) A line beginning with '# ----' which signals the end of the variable list

   4) Any number of lines containing data points, in the format:

      'Name_of_point var1 var2 var3 .... varN'

      all variable values should be integers

   Throughout, it is assumed that varN is the result

 - *.dat files* contain the same information as .qdat files, but the variable
   values can be anything (floats, ints, strings).  **These files should
   still contain quant_bounds!**

 - *.qdat.pkl file* contain a pickled (binary) representation of
   the data read in.  They stores, in order:

    1) A python list of the variable names

    2) A python list of lists with the quantization bounds

    3) A python list of the point names

    4) A python list of lists with the data points

"""
from __future__ import annotations
import rdkit.ML.Data.DataUtils
import typing
import csv
import numpy
import pickle
import rdkit.DataStructs.BitUtils
import rdkit.ML.Data.MLData
import rdkit.RDRandom
import rdkit.utils.fileutils
import re
_Shape = typing.Tuple[int, ...]

__all__ = [
    "BitUtils",
    "BuildDataSet",
    "BuildQuantDataSet",
    "CalcNPossibleUsingMap",
    "CountResults",
    "DBToData",
    "FilterData",
    "InitRandomNumbers",
    "MLData",
    "RandomizeActivities",
    "ReadGeneralExamples",
    "ReadQuantExamples",
    "ReadVars",
    "TakeEnsemble",
    "TextFileToData",
    "TextToData",
    "WriteData",
    "WritePickledData",
    "csv",
    "fileutils",
    "numpy",
    "permutation",
    "pickle",
    "random",
    "re"
]


