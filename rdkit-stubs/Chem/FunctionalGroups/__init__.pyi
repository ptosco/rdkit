from __future__ import annotations
import rdkit.Chem.FunctionalGroups
import typing
from _io import StringIO
import os
import rdkit.Chem
import rdkit.RDConfig
import re
import weakref

__all__ = [
    "BuildFuncGroupHierarchy",
    "Chem",
    "CreateMolFingerprint",
    "FGHierarchyNode",
    "FuncGroupFileParseError",
    "RDConfig",
    "StringIO",
    "groupDefns",
    "hierarchy",
    "lastData",
    "lastFilename",
    "os",
    "re",
    "weakref"
]


class FGHierarchyNode():
    children = None
    label = ''
    name = ''
    parent = None
    pattern = None
    removalReaction = None
    rxnSmarts = ''
    smarts = ''
    pass
class FuncGroupFileParseError(ValueError, Exception, BaseException):
    pass
groupDefns = {}
hierarchy = None
lastData = None
lastFilename = None
