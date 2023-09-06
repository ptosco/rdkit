""" functionality for generating an image showing the results of a composite model
voting on a data set

  Uses *Numeric* and *PIL*

"""
from __future__ import annotations
import rdkit.ML.MLUtils.VoteImg
import typing
import PIL.Image
import PIL.ImageDraw
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "BuildVoteImage",
    "CollectVotes",
    "Image",
    "ImageDraw",
    "Usage",
    "VoteAndBuildImage",
    "numpy"
]


