""" Defines the class _DecTreeNode_, used to represent decision trees

  _DecTreeNode_ is derived from _Tree.TreeNode_

"""
from __future__ import annotations
import rdkit.ML.DecTree.DecTree
import typing
import rdkit.ML.DecTree.Tree

__all__ = [
    "DecTreeNode",
    "Tree"
]


class DecTreeNode(rdkit.ML.DecTree.Tree.TreeNode):
    """
    This is used to represent decision trees

      _DecTreeNode_s are simultaneously the roots and branches of decision trees.
      Everything is nice and recursive.

      _DecTreeNode_s can save the following pieces of internal state, accessible via
        standard setter/getter functions:

        1) _Examples_: a list of examples which have been classified

        2) _BadExamples_: a list of examples which have been misclassified

        3) _TrainingExamples_: the list of examples used to train the tree

        4) _TestExamples_: the list of examples used to test the tree

     
    """
    pass
