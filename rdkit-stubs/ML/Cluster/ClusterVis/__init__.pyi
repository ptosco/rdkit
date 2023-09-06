"""Cluster tree visualization using Sping

"""
from __future__ import annotations
import rdkit.ML.Cluster.ClusterVis
import typing
import numpy
import rdkit.ML.Cluster.ClusterUtils
import rdkit.sping.colors
import rdkit.sping.pid
_Shape = typing.Tuple[int, ...]

__all__ = [
    "ClusterRenderer",
    "ClusterToImg",
    "ClusterToPDF",
    "ClusterToSVG",
    "ClusterUtils",
    "DrawClusterTree",
    "VisOpts",
    "numpy",
    "pid",
    "piddle"
]


class ClusterRenderer():
    pass
class VisOpts():
    """
    stores visualization options for cluster viewing

       **Instance variables**

         - x/yOffset: amount by which the drawing is offset from the edges of the canvas

         - lineColor: default color for drawing the cluster tree

         - lineWidth: the width of the lines used to draw the tree

     
    """
    hideColor: rdkit.sping.colors.Color # value = Color(0.80,0.80,0.80)
    hideWidth = 1.1
    highlightColor: rdkit.sping.colors.Color # value = Color(1.00,1.00,0.40)
    highlightRad = 10
    lineColor: rdkit.sping.colors.Color # value = Color(0.00,0.00,0.00)
    lineWidth = 2
    nodeColor: rdkit.sping.colors.Color # value = Color(1.00,0.40,0.40)
    nodeRad = 15
    terminalColors: list # value = [Color(1.00,0.00,0.00), Color(0.00,0.00,1.00), Color(1.00,1.00,0.00), Color(0.00,0.50,0.50), Color(0.00,0.80,0.00), Color(0.50,0.50,0.50), Color(0.80,0.30,0.30), Color(0.30,0.30,0.80), Color(0.80,0.80,0.30), Color(0.30,0.80,0.80)]
    xOffset = 20
    yOffset = 20
    pass
