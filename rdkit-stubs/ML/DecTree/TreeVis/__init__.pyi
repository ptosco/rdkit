""" functionality for drawing trees on sping canvases

"""
from __future__ import annotations
import rdkit.ML.DecTree.TreeVis
import typing
import math
import rdkit.sping.colors
import rdkit.sping.pid

__all__ = [
    "CalcTreeNodeSizes",
    "CalcTreeWidth",
    "DrawTree",
    "DrawTreeNode",
    "ResetTree",
    "SetNodeScales",
    "VisOpts",
    "math",
    "piddle",
    "visOpts"
]


class VisOpts():
    circColor: rdkit.sping.colors.Color # value = Color(0.60,0.60,0.90)
    circRad = 10
    highlightColor: rdkit.sping.colors.Color # value = Color(1.00,1.00,0.40)
    highlightWidth = 2
    horizOffset = 10
    labelFont: rdkit.sping.pid.Font # value = Font(10,0,0,0,'helvetica')
    lineColor: rdkit.sping.colors.Color # value = Color(0.00,0.00,0.00)
    lineWidth = 2
    maxCircRad = 16
    minCircRad = 4
    outlineColor: rdkit.sping.colors.Color # value = Color(-1.00,-1.00,-1.00)
    terminalEmptyColor: rdkit.sping.colors.Color # value = Color(0.80,0.80,0.20)
    terminalOffColor: rdkit.sping.colors.Color # value = Color(0.20,0.20,0.20)
    terminalOnColor: rdkit.sping.colors.Color # value = Color(0.80,0.80,0.80)
    vertOffset = 50
    pass
visOpts: rdkit.ML.DecTree.TreeVis.VisOpts
