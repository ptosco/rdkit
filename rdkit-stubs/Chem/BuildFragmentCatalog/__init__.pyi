"""  command line utility for working with FragmentCatalogs (CASE-type analysis)

**Usage**

  BuildFragmentCatalog [optional args] <filename>

 filename, the name of a delimited text file containing InData, is required
 for some modes of operation (see below)

**Command Line Arguments**

 - -n *maxNumMols*:  specify the maximum number of molecules to be processed

 - -b: build the catalog and OnBitLists
    *requires InData*

 - -s: score compounds
    *requires InData and a Catalog, can use OnBitLists*

 - -g: calculate info gains
    *requires Scores*

 - -d: show details about high-ranking fragments
    *requires a Catalog and Gains*

 - --catalog=*filename*: filename with the pickled catalog.
    If -b is provided, this file will be overwritten.

 - --onbits=*filename*: filename to hold the pickled OnBitLists.
   If -b is provided, this file will be overwritten

 - --scores=*filename*: filename to hold the text score data.
   If -s is provided, this file will be overwritten

 - --gains=*filename*: filename to hold the text gains data.
   If -g is provided, this file will be overwritten

 - --details=*filename*: filename to hold the text details data.
   If -d is provided, this file will be overwritten.

 - --minPath=2: specify the minimum length for a path

 - --maxPath=6: specify the maximum length for a path

 - --smiCol=1: specify which column in the input data file contains
     SMILES

 - --actCol=-1: specify which column in the input data file contains
     activities

 - --nActs=2: specify the number of possible activity values

 - --nBits=-1: specify the maximum number of bits to show details for

"""
from __future__ import annotations
import rdkit.Chem.BuildFragmentCatalog
import typing
from rdkit.Dbase.DbConnection import DbConnect
import numpy
import os
import pickle
import rdkit.Chem.FragmentCatalog
import rdkit.ML.InfoTheory
import rdkit.RDConfig
import sys
_Shape = typing.Tuple[int, ...]

__all__ = [
    "BuildCatalog",
    "CalcGains",
    "CalcGainsFromFps",
    "DbConnect",
    "FragmentCatalog",
    "InfoTheory",
    "OutputGainsData",
    "ParseArgs",
    "ProcessGainsData",
    "RDConfig",
    "RunDetails",
    "ScoreFromLists",
    "ScoreMolecules",
    "ShowDetails",
    "SupplierFromDetails",
    "Usage",
    "message",
    "numpy",
    "os",
    "pickle",
    "sys"
]


class RunDetails():
    actCol = -1
    biasList = None
    catalogName = None
    dbName = ''
    delim = ','
    detailsName = None
    doBuild = 0
    doDetails = 0
    doGains = 0
    doScore = 0
    doSigs = 0
    fpName = None
    gainsName = None
    hasTitle = 1
    inFileName = None
    maxPath = 6
    minPath = 2
    nActs = 2
    nBits = -1
    nameCol = -1
    numMols = -1
    onBitsName = None
    scoresName = None
    smiCol = 1
    tableName = None
    topN = -1
    pass
