""" Define the class _KNNClassificationModel_, used to represent a k-nearest neighbhors
    classification model

    Inherits from _KNNModel_
"""
from __future__ import annotations
import rdkit.ML.KNN.KNNClassificationModel
import typing
import rdkit.ML.KNN.KNNModel

__all__ = [
    "KNNClassificationModel",
    "KNNModel"
]


class KNNClassificationModel(rdkit.ML.KNN.KNNModel.KNNModel):
    """
    This is used to represent a k-nearest neighbor classifier

     
    """
    pass
