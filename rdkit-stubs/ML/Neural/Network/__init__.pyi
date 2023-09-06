""" Contains the class _Network_ which is used to represent neural nets

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
import rdkit.ML.Neural.Network
import typing
import numpy
import random
import rdkit.ML.Neural.ActFuncs
import rdkit.ML.Neural.NetNode
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ActFuncs",
    "NetNode",
    "Network",
    "numpy",
    "random"
]


class Network():
    """
    a neural network

     
    """
    pass
