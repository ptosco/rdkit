""" Defines the class _QuantTreeNode_, used to represent decision trees with automatic
 quantization bounds

  _QuantTreeNode_ is derived from _DecTree.DecTreeNode_

"""
from __future__ import annotations
import rdkit.ML.DecTree.QuantTree
import typing
import rdkit.ML.DecTree.DecTree
import rdkit.ML.DecTree.Tree

__all__ = [
    "DecTree",
    "QuantTreeNode",
    "Tree"
]


class QuantTreeNode(rdkit.ML.DecTree.DecTree.DecTreeNode, rdkit.ML.DecTree.Tree.TreeNode):
    """
      
    """
    __hash__ = None
    pass
