""" Training algorithms for feed-forward neural nets

  Unless noted otherwise, algorithms and notation are taken from:
  "Artificial Neural Networks: Theory and Applications",
    Dan W. Patterson, Prentice Hall, 1996

"""
from __future__ import annotations
import rdkit.ML.Neural.Trainers
import typing
import numpy
_Shape = typing.Tuple[int, ...]

__all__ = [
    "BackProp",
    "Trainer",
    "numpy"
]


class Trainer():
    """
    "virtual base class" for network trainers

     
    """
    pass
class BackProp(Trainer):
    """
    implement back propagation (algorithm on pp 153-154 of Patterson)

       I don't *think* that I've made any assumptions about the connectivity of
         the net (i.e. full connectivity between layers is not required).

       **NOTE:** this code is currently making the assumption that the activation
         functions on the nodes in the network are capable of calculating their
         derivatives using only their values (i.e. a DerivFromVal method should
         exist).  This shouldn't be too hard to change.

      
    """
    pass
