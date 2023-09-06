from __future__ import annotations
import rdkit.VLib.NodeLib.SDSupply
import typing
from rdkit.VLib.Supply import SupplyNode
import rdkit.Chem
import rdkit.VLib.Node
import rdkit.VLib.Supply

__all__ = [
    "Chem",
    "SDSupplyNode",
    "SupplyNode"
]


class SDSupplyNode(rdkit.VLib.Supply.SupplyNode, rdkit.VLib.Node.VLibNode):
    """
    SD supplier

       Sample Usage:
         >>> import os
         >>> from rdkit import RDConfig
         >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',                               'test_data','NCI_aids.10.sdf')
         >>> suppl = SDSupplyNode(fileN)
         >>> ms = [x for x in suppl]
         >>> len(ms)
         10
         >>> ms[0].GetProp("_Name")
         '48'
         >>> ms[1].GetProp("_Name")
         '78'
         >>> suppl.reset()
         >>> suppl.next().GetProp("_Name")
         '48'
         >>> suppl.next().GetProp("_Name")
         '78'

       
    """
    pass
