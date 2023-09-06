""" Generic file manipulation stuff

"""
from __future__ import annotations
import rdkit.ML.files
import typing
import numpy
import re
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ReFile",
    "ReadDataFile",
    "numpy",
    "re"
]


class ReFile():
    """
    convenience class for dealing with files with comments

      blank (all whitespace) lines, and lines beginning with comment
        characters are skipped.

      anything following a comment character on a line is stripped off
      
    """
    pass
