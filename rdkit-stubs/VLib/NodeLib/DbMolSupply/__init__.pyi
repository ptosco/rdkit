from __future__ import annotations
import rdkit.VLib.NodeLib.DbMolSupply
import typing
from rdkit.VLib.Supply import SupplyNode
import os
import rdkit.Chem
import rdkit.Chem.Suppliers.DbMolSupplier
import rdkit.RDConfig
import rdkit.VLib.Node
import rdkit.VLib.Supply
import sys

__all__ = [
    "Chem",
    "DbMolSupplier",
    "DbMolSupplyNode",
    "GetNode",
    "RDConfig",
    "SupplyNode",
    "os",
    "sys"
]


class DbMolSupplyNode(rdkit.VLib.Supply.SupplyNode, rdkit.VLib.Node.VLibNode):
    """
    Supplies molecules from a db result set:

     Sample Usage:
       >>> from rdkit.Dbase.DbConnection import DbConnect
       >>> dbName = os.path.join(RDConfig.RDCodeDir,'Chem','Fingerprints',                             'test_data','data.gdb')
       >>> conn = DbConnect(dbName,'simple_mols')
       >>> dataset = conn.GetData()
       >>> suppl = DbMolSupplyNode(dataset)
       >>> ms = [x for x in suppl]
       >>> len(ms)
       12
       >>> ms[0].GetProp("ID")
       'ether-1'
       >>> ms[10].GetProp("ID")
       'acid-4'
       >>> suppl.reset()
       >>> suppl.next().GetProp("ID")
       'ether-1'
       >>> suppl.next().GetProp("ID")
       'acid-1'
       >>> suppl.reset()
     
     
    """
    pass
