""" Defines the class SigTreeNode, used to represent trees that
 use signatures (bit vectors) to represent data.  As inputs (examples),
 SigTreeNode's expect 3-sequences: (label,sig,act)

  _SigTreeNode_ is derived from _DecTree.DecTreeNode_

"""
from __future__ import annotations
import rdkit.ML.DecTree.SigTree
import typing
from rdkit.DataStructs.VectCollection import VectCollection
import copy
import rdkit.ML.DecTree.DecTree
import rdkit.ML.DecTree.Tree

__all__ = [
    "DecTree",
    "SigTreeNode",
    "VectCollection",
    "copy"
]


class SigTreeNode(rdkit.ML.DecTree.DecTree.DecTreeNode, rdkit.ML.DecTree.Tree.TreeNode):
    """
      
    """
    pass
