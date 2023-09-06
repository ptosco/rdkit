""" Configuration for the RDKit Python code

"""
from __future__ import annotations
import rdkit.RDConfig
import typing
import os
import sqlite3
import sys

__all__ = [
    "ObsoleteCodeError",
    "RDBaseDir",
    "RDBinDir",
    "RDCodeDir",
    "RDContribDir",
    "RDDataDatabase",
    "RDDataDir",
    "RDDemoDir",
    "RDDocsDir",
    "RDProjDir",
    "RDTestDatabase",
    "UnimplementedCodeError",
    "defaultDBPassword",
    "defaultDBUser",
    "molViewer",
    "os",
    "pythonExe",
    "pythonTestCommand",
    "rpcTestPort",
    "sqlite3",
    "sys",
    "usePgSQL",
    "useSqlLite"
]


class ObsoleteCodeError(Exception, BaseException):
    pass
class UnimplementedCodeError(Exception, BaseException):
    pass
RDBaseDir = '/scratch/toscopa1/src/rdkit'
RDBinDir = '/scratch/toscopa1/src/rdkit/bin'
RDCodeDir = '/scratch/toscopa1/src/rdkit/rdkit'
RDContribDir = '/scratch/toscopa1/src/rdkit/Contrib'
RDDataDatabase = '/scratch/toscopa1/src/rdkit/Data/RDData.sqlt'
RDDataDir = '/scratch/toscopa1/src/rdkit/Data'
RDDemoDir = '/scratch/toscopa1/src/rdkit/Demo'
RDDocsDir = '/scratch/toscopa1/src/rdkit/Docs'
RDProjDir = '/scratch/toscopa1/src/rdkit/Projects'
RDTestDatabase = '/scratch/toscopa1/src/rdkit/Data/RDTests.sqlt'
defaultDBPassword = 'masterkey'
defaultDBUser = 'sysdba'
molViewer = 'PYMOL'
pythonExe = '/scratch/toscopa1/.conda/envs/pyds_v1.2L/bin/python'
pythonTestCommand = 'python'
rpcTestPort = 8423
usePgSQL = False
useSqlLite = True
