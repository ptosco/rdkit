""" command line utility to report on the contributions of descriptors to
tree-based composite models

Usage:  AnalyzeComposite [optional args] <models>

      <models>: file name(s) of pickled composite model(s)
        (this is the name of the db table if using a database)

    Optional Arguments:

      -n number: the number of levels of each model to consider

      -d dbname: the database from which to read the models

      -N Note: the note string to search for to pull models from the database

      -v: be verbose whilst screening
"""
from __future__ import annotations
import rdkit.ML.AnalyzeComposite
import typing
from rdkit.Dbase.DbConnection import DbConnect
import numpy
import pickle
import rdkit.ML.Data.Stats
import rdkit.ML.DecTree.Tree
import rdkit.ML.DecTree.TreeUtils
import rdkit.ML.ScreenComposite
import sys
_Shape = typing.Tuple[int, ...]

__all__ = [
    "DbConnect",
    "ErrorStats",
    "ProcessIt",
    "ScreenComposite",
    "ShowStats",
    "Stats",
    "Tree",
    "TreeUtils",
    "Usage",
    "numpy",
    "pickle",
    "sys"
]


__VERSION_STRING = '2.2.0'
