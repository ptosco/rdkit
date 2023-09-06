from __future__ import annotations
import rdkit.VLib.NodeLib.SmilesOutput
import typing
import rdkit.Chem
import rdkit.VLib.Node

__all__ = [
    "BaseOutputNode",
    "Chem",
    "OutputNode"
]


class OutputNode(rdkit.VLib.Node.VLibNode):
    """
    dumps smiles output

       Assumptions:

         - destination supports a write() method

         - inputs (parents) can be stepped through in lockstep


       Usage Example:
         >>> smis = ['C1CCC1','C1CC1','C=O','NCC']
         >>> mols = [Chem.MolFromSmiles(x) for x in smis]
         >>> from rdkit.VLib.Supply import SupplyNode
         >>> suppl = SupplyNode(contents=mols)
         >>> from io import StringIO
         >>> sio = StringIO()
         >>> node = OutputNode(dest=sio,delim=', ')
         >>> node.AddParent(suppl)
         >>> ms = [x for x in node]
         >>> len(ms)
         4
         >>> txt = sio.getvalue()
         >>> repr(txt)
         "'1, C1CCC1\\n2, C1CC1\\n3, C=O\\n4, CCN\\n'"

       
    """
    pass
BaseOutputNode = rdkit.VLib.Output.OutputNode
