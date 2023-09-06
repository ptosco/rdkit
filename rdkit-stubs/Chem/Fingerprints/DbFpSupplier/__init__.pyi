""" Supplies a class for working with fingerprints from databases
#DOC

"""
from __future__ import annotations
import rdkit.Chem.Fingerprints.DbFpSupplier
import typing
from rdkit.VLib.Node import VLibNode
import pickle
import rdkit.DataStructs
import rdkit.VLib.Node

__all__ = [
    "DataStructs",
    "DbFpSupplier",
    "ForwardDbFpSupplier",
    "RandomAccessDbFpSupplier",
    "VLibNode",
    "pickle"
]


class DbFpSupplier(rdkit.VLib.Node.VLibNode):
    """
    new fps come back with all additional fields from the
    database set in a "_fieldsFromDb" data member
    """
    pass
class ForwardDbFpSupplier(DbFpSupplier, rdkit.VLib.Node.VLibNode):
    """
    DbFp supplier supporting only forward iteration

       >>> from rdkit import RDConfig
       >>> from rdkit.Dbase.DbConnection import DbConnect
       >>> fName = RDConfig.RDTestDatabase
       >>> conn = DbConnect(fName,'simple_combined')
       >>> suppl = ForwardDbFpSupplier(conn.GetData())

       we can loop over the supplied fingerprints:
       
       >>> fps = []
       >>> for fp in suppl:
       ...   fps.append(fp)
       >>> len(fps)
       12

       
    """
    pass
class RandomAccessDbFpSupplier(DbFpSupplier, rdkit.VLib.Node.VLibNode):
    """
    DbFp supplier supporting random access:

     >>> import os.path
     >>> from rdkit import RDConfig
     >>> from rdkit.Dbase.DbConnection import DbConnect
     >>> fName = RDConfig.RDTestDatabase
     >>> conn = DbConnect(fName,'simple_combined')
     >>> suppl = RandomAccessDbFpSupplier(conn.GetData())
     >>> len(suppl)
     12

     we can pull individual fingerprints:

     >>> fp = suppl[5]
     >>> fp.GetNumBits()
     128
     >>> fp.GetNumOnBits()
     54

     a standard loop over the fingerprints:

     >>> fps = []
     >>> for fp in suppl:
     ...   fps.append(fp)
     >>> len(fps)
     12

     or we can use an indexed loop:

     >>> fps = [None] * len(suppl)
     >>> for i in range(len(suppl)):
     ...   fps[i] = suppl[i]
     >>> len(fps)
     12

     
    """
    pass
