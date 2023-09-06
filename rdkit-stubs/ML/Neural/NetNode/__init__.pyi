""" Contains the class _NetNode_ which is used to represent nodes in neural nets

**Network Architecture:**

  A tacit assumption in all of this stuff is that we're dealing with
  feedforward networks.

  The network itself is stored as a list of _NetNode_ objects.  The list
  is ordered in the sense that nodes in earlier/later layers than a
  given node are guaranteed to come before/after that node in the list.
  This way we can easily generate the values of each node by moving
  sequentially through the list, we're guaranteed that every input for a
  node has already been filled in.

  Each node stores a list (_inputNodes_) of indices of its inputs in the
  main node list.

"""
from __future__ import annotations
import rdkit.ML.Neural.NetNode
import typing
import numpy
import rdkit.ML.Neural.ActFuncs
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ActFuncs",
    "NetNode",
    "numpy"
]


class NetNode():
    """
    a node in a neural network

     
    """
    pass
