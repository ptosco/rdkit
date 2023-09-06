""" uses pymol to interact with molecules

"""
from __future__ import annotations
import rdkit.Chem.PyMol
import typing
import os
import rdkit.Chem
import sys
import tempfile

__all__ = [
    "Chem",
    "MolViewer",
    "Server",
    "os",
    "sys",
    "tempfile"
]


class MolViewer():
    pass
_server = None
Server = xmlrpc.client.ServerProxy
