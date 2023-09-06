""" Various storage (molecular and otherwise) functionality

"""
from __future__ import annotations
import rdkit.Dbase.StorageUtils
import typing
import rdkit.Dbase.DbModule
import rdkit.RDConfig

__all__ = [
    "DbModule",
    "GetNextId",
    "GetNextRDId",
    "IndexToRDId",
    "RDConfig",
    "RDIdToInt",
    "RegisterItem",
    "RegisterItems",
    "ValidateRDId"
]


__test__ = {'roundtrip': "\n>>> ValidateRDId(IndexToRDId(100))\n1\n>>> ValidateRDId(IndexToRDId(10000,leadText='foo'))\n1\n>>> indices = [1,100,1000,1000000]\n>>> vals = []\n>>> for idx in indices:\n...   vals.append(RDIdToInt(IndexToRDId(idx)))\n>>> vals == indices\n1\n\n"}
_roundtripTests = "\n>>> ValidateRDId(IndexToRDId(100))\n1\n>>> ValidateRDId(IndexToRDId(10000,leadText='foo'))\n1\n>>> indices = [1,100,1000,1000000]\n>>> vals = []\n>>> for idx in indices:\n...   vals.append(RDIdToInt(IndexToRDId(idx)))\n>>> vals == indices\n1\n\n"
