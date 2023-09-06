from __future__ import annotations
import rdkit.VLib.NodeLib.SmilesDupeFilter
import typing
from rdkit.VLib.Filter import FilterNode
import rdkit.Chem
import rdkit.VLib.Filter
import rdkit.VLib.Node

__all__ = [
    "Chem",
    "DupeFilter",
    "FilterNode"
]


class DupeFilter(rdkit.VLib.Filter.FilterNode, rdkit.VLib.Node.VLibNode):
    """
    canonical-smiles based duplicate filter

     Assumptions:

       - inputs are molecules


     Sample Usage:
       >>> import os
       >>> from rdkit import RDConfig
       >>> from rdkit.VLib.NodeLib.SDSupply import SDSupplyNode
       >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',                             'test_data','NCI_aids.10.sdf')
       >>> suppl = SDSupplyNode(fileN)
       >>> filt = DupeFilter()
       >>> filt.AddParent(suppl)
       >>> ms = [x for x in filt]
       >>> len(ms)
       10
       >>> ms[0].GetProp("_Name")
       '48'
       >>> ms[1].GetProp("_Name")
       '78'
       >>> filt.reset()
       >>> filt.next().GetProp("_Name")
       '48'


     
    """
    pass
