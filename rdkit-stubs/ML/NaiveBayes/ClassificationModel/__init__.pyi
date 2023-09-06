""" Defines Naive Baysean classification model
   Based on development in: Chapter 6 of "Machine Learning" by Tom Mitchell

"""
from __future__ import annotations
import rdkit.ML.NaiveBayes.ClassificationModel
import typing
import numpy
import rdkit.ML.Data.Quantize
_Shape = typing.Tuple[int, ...]

__all__ = [
    "NaiveBayesClassifier",
    "Quantize",
    "numpy"
]


class NaiveBayesClassifier():
    """
    _NaiveBayesClassifier_s can save the following pieces of internal state, accessible via
    standard setter/getter functions:

    1) _Examples_: a list of examples which have been predicted

    2) _TrainingExamples_: List of training examples - the descriptor value of these examples
      are quantized based on info gain using ML/Data/Quantize.py if necessary

    3) _TestExamples_: the list of examples used to test the model

    4) _BadExamples_ : list of examples that were incorrectly classified

    4) _QBoundVals_: Quant bound values for each varaible - a list of lists

    5) _QBounds_ : Number of bounds for each variable
    """
    pass
