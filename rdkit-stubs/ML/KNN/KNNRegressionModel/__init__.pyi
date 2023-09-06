""" Define the class _KNNRegressionModel_, used to represent a k-nearest neighbhors
regression model

    Inherits from _KNNModel_
"""
from __future__ import annotations
import rdkit.ML.KNN.KNNRegressionModel
import typing
import rdkit.ML.KNN.KNNModel

__all__ = [
    "KNNModel",
    "KNNRegressionModel"
]


class KNNRegressionModel(rdkit.ML.KNN.KNNModel.KNNModel):
    """
    This is used to represent a k-nearest neighbor classifier

     
    """
    pass
