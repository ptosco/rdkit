""" Activation functions for neural network nodes

Activation functions should implement the following API:

 - _Eval(x)_: returns the value of the function at a given point

 - _Deriv(x)_: returns the derivative of the function at a given point

The current Backprop implementation also requires:

 - _DerivFromVal(val)_: returns the derivative of the function when its
                        value is val

In all cases _x_ is a float as is the value returned.

"""
from __future__ import annotations
import rdkit.ML.Neural.ActFuncs
import typing
import math

__all__ = [
    "ActFunc",
    "Sigmoid",
    "TanH",
    "math"
]


class ActFunc():
    """
    "virtual base class" for activation functions

     
    """
    pass
class Sigmoid(ActFunc):
    """
    the standard sigmoidal function 
    """
    pass
class TanH(ActFunc):
    """
    the standard hyperbolic tangent function 
    """
    pass
