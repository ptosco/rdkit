from __future__ import annotations
import rdkit.ML.ModelPackage.Packager
import typing

__all__ = [
    "ClassificationError",
    "DescriptorCalculationError",
    "ModelPackage"
]


class ClassificationError(Exception, BaseException):
    """
    used to signal problems generating predictions 
    """
    pass
class DescriptorCalculationError(Exception, BaseException):
    """
    used to signal problems generating descriptor values 
    """
    pass
class ModelPackage():
    """
    a container class to package a composite model with a descriptor
     calculator so that objects needing predictions (compounds, molecules, etc.)
     can be passed directly in without worrying about generating descriptors

     
    """
    pass
